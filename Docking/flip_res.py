import os
import numpy as np
from modeller import *
from openmm import *
from openmm.app import *
from pdbfixer import PDBFixer
from scipy.spatial import cKDTree

try:
    from simtk import unit
except ImportError:
    from openmm import unit

# noinspection PyUnresolvedReferences
nanometer = unit.nanometer


# ResidueFlipper – wykonuje flipy ASN, GLN i HIS.
class ResidueFlipper:

    def __init__(self, input_pdb_path, aln=None):
        """
        Przyjmuje ścieżkę do PDB i alignment.
        Inicjalizuje środowisko Modeller.
        """
        self.input_pdb_path = input_pdb_path
        self.aln = aln
        self.env = environ()
        self.env.io.hetatm = True
        self.output_pdb_path = os.path.splitext(input_pdb_path)[0] + "_flipped_openmm.pdb"

    def count_local_hbonds(self, residue, topology, positions, cutoff=0.35):
        """
        Liczy lokalne wiązania wodorowe O/N wokół reszty.
        Wykorzystuje KDTree dla efektywności.
        """
        atoms = list(topology.atoms())
        pos_array = np.array([[positions[i].x, positions[i].y, positions[i].z] for i in range(len(atoms))])
        elements = [atom.element.symbol for atom in atoms]

        # tylko O i N
        donor_acceptor_indices = [i for i, el in enumerate(elements) if el in ('O', 'N')]
        donor_acceptor_coords = pos_array[donor_acceptor_indices]

        kdtree = cKDTree(donor_acceptor_coords)

        res_indices = [i for i, atom in enumerate(atoms) if atom.residue == residue and elements[i] in ('O', 'N')]
        res_coords = pos_array[res_indices]

        # Znajduje wszystkie punkty w zasięgu koordynatów
        count = 0
        for coord in res_coords:
            close = kdtree.query_ball_point(coord, cutoff)
            count += len(close)

        return count

    def flip_amide_single(self, residue, structure, positions, resname):
        """
        Flip ASN/GLN: obraca grupę amidową (O, N, hydrogens) o 180°
        wokół odpowiedniego wiązania (ASN: CG–CB, GLN: CD–CG).

        Modyfikuje obiekt 'positions' w miejscu.
        """
        try:
            atom_dict = {atom.name: atom for atom in residue.atoms()}

            if resname == 'ASN':
                required = ['CG', 'OD1', 'ND2']
                pivot_name = 'CG'
                axis_partner = 'CB'  # fallback: CA
                group_base = ['OD1', 'ND2']
                hyd_prefix = 'HD2'
            elif resname == 'GLN':
                required = ['CD', 'OE1', 'NE2']
                pivot_name = 'CD'
                axis_partner = 'CG'  # fallback: CB
                group_base = ['OE1', 'NE2']
                hyd_prefix = 'HE2'
            else:
                return False

            # Sprawdź wymagane atomy
            if not all(name in atom_dict for name in required):
                # Dzieje się tak, gdy PDBFixer doda brakujące atomy (np. w GLN 222),
                # ale plik wejściowy był tak zepsuty, że nie dodał ich poprawnie.
                print(f"Ostrzeżenie: Pomijanie flipa {resname} {residue.id} - brak wymaganych atomów ({required})")
                return False

            # Mapowanie na indeksy
            atom_lookup = {atom.name: i for i, atom in enumerate(structure.atoms()) if atom.residue == residue}

            # Pivot atom
            pivot_idx = atom_lookup[pivot_name]
            pivot_pos = positions[pivot_idx].value_in_unit(nanometer)
            center = np.array([pivot_pos.x, pivot_pos.y, pivot_pos.z])

            # Axis partner (np. CB lub CG)
            if axis_partner in atom_lookup:
                partner_idx = atom_lookup[axis_partner]
            else:
                partner_idx = atom_lookup['CA']  # fallback
            partner_pos = positions[partner_idx].value_in_unit(nanometer)
            axis = np.array([partner_pos.x, partner_pos.y, partner_pos.z]) - center

            # Gdyby CD
            if np.linalg.norm(axis) < 1e-6:
                print(f"Ostrzeżenie: Pomijanie flipa {resname} {residue.id} - oś rotacji ma zerową długość.")
                return False

            axis /= np.linalg.norm(axis)

            # Grupa do obrotu
            group_names = group_base + [a.name for a in residue.atoms() if a.name.startswith(hyd_prefix)]
            group_indices = [atom_lookup[n] for n in group_names if n in atom_lookup]

            group_positions = [
                np.array([positions[i].x, positions[i].y, positions[i].z])
                for i in group_indices
            ]

            # Obrót o 180°
            new_positions = self.rotate_atoms(group_positions, center, axis, np.pi)

            # Zapisz nowe współrzędne
            for i, idx in enumerate(group_indices):
                new_pos = new_positions[i]
                positions[idx] = Vec3(new_pos[0], new_pos[1], new_pos[2]) * nanometer

            return True
        except Exception as e:
            print(f"Błąd w flip_amide_single ({resname} {residue.id}): {e}")
            return False

    def flip_HIS_single(self, residue, structure, positions):
        """
        Flip HIS: obraca cały pierścień (i jego wodory) o 180°
        wokół osi CB-CG.
        Modyfikuje obiekt 'positions' w miejscu.
        """
        try:
            atom_dict = {atom.name: atom for atom in residue.atoms()}

            # Wymagamy atomów do zdefiniowania osi (CB, CG) i pierścienia
            required = ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']
            if not all(name in atom_dict for name in required):
                print(f"Ostrzeżenie: Pomijanie flipa HIS {residue.id} - brak wymaganych atomów ({required})")
                return False

            # Mapowanie na indeksy LOKALNE DLA TEJ RESZTY
            atom_lookup = {atom.name: i for i, atom in enumerate(structure.atoms()) if atom.residue == residue}

            # Definicja osi rotacji: CB -> CG
            cb_idx = atom_lookup['CB']
            cg_idx = atom_lookup['CG']

            cb_pos = positions[cb_idx].value_in_unit(nanometer)
            cg_pos = positions[cg_idx].value_in_unit(nanometer)

            center = cb_pos # Pivot (punkt na osi) to CB
            axis = cg_pos - center # Wektor osi

            if np.linalg.norm(axis) < 1e-6:
                print(f"Ostrzeżenie: Pomijanie flipa HIS {residue.id} - oś rotacji CB-CG ma zerową długość.")
                return False

            axis /= np.linalg.norm(axis)

            # Grupa do rotacji: Cały pierścień (włącznie z CG)
            ring_atom_names = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']

            # Dodaj wodory pierścienia do rotacji
            hyd_names = [a.name for a in residue.atoms() if a.name.startswith('H') and a.name not in ('H', 'HA', 'HB2', 'HB3')]
            ring_atom_names.extend(hyd_names)

            ring_indices = [atom_lookup[name] for name in ring_atom_names if name in atom_lookup]

            ring_positions = [
                np.array([positions[i].x, positions[i].y, positions[i].z])
                for i in ring_indices
            ]

            # Obrót o 180°
            new_positions = self.rotate_atoms(ring_positions, center, axis, np.pi)

            for i, idx in enumerate(ring_indices):
                new_pos = new_positions[i]
                positions[idx] = Vec3(new_pos[0], new_pos[1], new_pos[2]) * nanometer
            return True
        except Exception as e:
            print(f"Błąd w flip_HIS_single (HIS {residue.id}): {e}")
            return False

    def rotate_atoms(self, atom_positions, center, axis, angle):
        """
        Właściwa funkcja do obrotu pozycji atomów.
        """
        axis_norm = np.linalg.norm(axis)
        if axis_norm < 1e-6:
            return atom_positions # Nie można obracać wokół osi zerowej

        axis = axis / axis_norm
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        rotated_positions = []
        for pos in atom_positions:
            rel_pos = pos - center
            rotated = (rel_pos * cos_theta +
                       np.cross(axis, rel_pos) * sin_theta +
                       axis * np.dot(axis, rel_pos) * (1 - cos_theta))
            rotated_positions.append(rotated + center)

        return rotated_positions

    def run(self):
        """
        Główna funkcja:
            - wczytuje strukturę,
            - testuje flipy na KOPII współrzędnych,
            - zatwierdza flip tylko jeśli jest zysk HBonds,
            - zapisuje wynik.
        """
        # Wczytanie pliku
        fixer = PDBFixer(filename=self.input_pdb_path)
        structure = fixer.topology
        positions = fixer.positions

        flip_count = 0

        # Przechodzimy po resztach i flipujemy jeśli jest to potrzebne
        for chain in structure.chains():
            for residue in chain.residues():
                resname = residue.name.upper()
                resid = residue.id

                if resname not in ['ASN', 'GLN', 'HIS']:
                    continue

                # 1. Policz wiązania na ORYGINALNYCH pozycjach
                hbonds_before = self.count_local_hbonds(residue, structure, positions)

                # 2. Stwórz tymczasową KOPIĘ do testowania flipa
                temp_positions = positions[:]
                flipped = False

                # 3. Wykonaj flip na KOPII
                if resname in ['ASN', 'GLN']:
                    flipped = self.flip_amide_single(residue, structure, temp_positions, resname)
                elif resname == 'HIS':
                    flipped = self.flip_HIS_single(residue, structure, temp_positions)

                # Jeśli doszło do flipa, rotacja atomów już została wykonana
                if not flipped:
                    continue

                # 4. Policz wiązania na sflipowanej KOPII
                hbonds_after = self.count_local_hbonds(residue, structure, temp_positions)

                # 5. Porównaj i podejmij decyzję
                if hbonds_after > hbonds_before:
                    # ZYSK: Nadpisz oryginalne pozycje sflipowaną kopią
                    print(f"Zastosowano flip {resname} {resid}, zwiększono lokalne wiązania wodorowe: {hbonds_before} → {hbonds_after}")
                    positions[:] = temp_positions[:] # Zastosuj zmiany do oryginału
                    flip_count += 1
                else:
                    # BRAK ZYSKU: Nie rób nic, kopia zostanie porzucona.
                    print(f"Nie zastosowano flipa {resname} {resid}, brak poprawy: {hbonds_before} → {hbonds_after}")

        # Zapis do pliku
        with open(self.output_pdb_path, 'w') as out_file:
            PDBFile.writeFile(structure, positions, out_file)

        print(f"\nZastosowano łącznie {flip_count} flipów.")
        print(f">>> Gotowe flipy zapisano do: {self.output_pdb_path}")

        return self.output_pdb_path

# Użycie:
# flipper = ResidueFlipper("sciezka/do/pliku.pdb")
# flipper.run()