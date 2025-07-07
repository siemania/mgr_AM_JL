import os
import numpy as np
from modeller import *
from modeller.automodel import *
from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
from scipy.spatial import cKDTree
try:
    from simtk import unit
except ImportError:
    from openmm import unit

nanometer = unit.nanometer

# ResidueFlipper – wykonuje flipy ASN, GLN i HIS.
class ResidueFlipper:
    def __init__(self, input_pdb_path, aln):
        """
        Przyjmuje ścieżkę do PDB i alignment.
        Inicjalizuje środowisko Modeller.
        """
        self.input_pdb_path = input_pdb_path
        self.aln = aln
        self.env = environ()
        self.env.io.hetatm = True

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

        count = 0
        for coord in res_coords:
            close = kdtree.query_ball_point(coord, cutoff)
            count += len(close)

        return count

    def flip_and_count_hbonds(self, resname, residue, structure, positions):
        """
        Próbuje flipować i sprawdza, czy liczba HBonds wzrosła.
        """
        original_positions = positions[:]
        flipped = False

        if resname == 'ASN':
            flipped = self.flip_ASN_single(residue, structure, positions)
        elif resname == 'GLN':
            flipped = self.flip_GLN_single(residue, structure, positions)
        elif resname == 'HIS':
            flipped = self.flip_HIS_single(residue, None, None, structure, positions)

        if not flipped:
            return original_positions, False

        return positions[:], True

    def flip_GLN_single(self, residue, structure, positions):
        """
        Flip GLN: zamienia OE1 i NE2.
        """
        try:
            atom_indices = {atom.name: atom.index for atom in structure.atoms() if atom.residue == residue}
            if 'OE1' not in atom_indices or 'NE2' not in atom_indices:
                return False
            oe1_idx = atom_indices['OE1']
            ne2_idx = atom_indices['NE2']
            positions[oe1_idx], positions[ne2_idx] = positions[ne2_idx], positions[oe1_idx]
            return True
        except:
            return False

    def flip_ASN_single(self, residue, structure, positions):
        """
        Flip ASN: zamienia OD1 i ND2.
        """
        try:
            atom_indices = {atom.name: atom.index for atom in structure.atoms() if atom.residue == residue}
            if 'OD1' not in atom_indices or 'ND2' not in atom_indices:
                return False
            od1_idx = atom_indices['OD1']
            nd2_idx = atom_indices['ND2']
            positions[od1_idx], positions[nd2_idx] = positions[nd2_idx], positions[od1_idx]
            return True
        except:
            return False

    def flip_HIS_single(self, residue, mdl, fixer, structure, positions):
        """
        Flip HIS: obraca pierścień o 180°.
        """
        try:
            atom_dict = {atom.name: atom for atom in residue.atoms()}
            if not all(name in atom_dict for name in ['CG', 'CD2', 'ND1', 'CE1', 'NE2']):
                return False
            atom_lookup = {atom.name: i for i, atom in enumerate(structure.atoms())}
            cg_idx = atom_lookup['CG']
            cd2_idx = atom_lookup['CD2']
            cg_pos = positions[cg_idx].value_in_unit(nanometer)
            cd2_pos = positions[cd2_idx].value_in_unit(nanometer)
            center = np.array([cg_pos.x, cg_pos.y, cg_pos.z])
            axis = np.array([cd2_pos.x, cd2_pos.y, cd2_pos.z]) - center
            axis /= np.linalg.norm(axis)
            ring_atom_names = ['ND1', 'CD2', 'CE1', 'NE2']
            ring_indices = [atom_lookup[name] for name in ring_atom_names]
            ring_positions = [
                np.array([positions[i].x, positions[i].y, positions[i].z])
                for i in ring_indices
            ]
            new_positions = self.rotate_atoms(ring_positions, center, axis, np.pi)
            for i, idx in enumerate(ring_indices):
                new_pos = new_positions[i]
                positions[idx] = Vec3(new_pos[0], new_pos[1], new_pos[2]) * nanometer
            return True
        except:
            return False

    def rotate_atoms(self, atom_positions, center, axis, angle):
        """
        Właściwa funkcja do obrotu pozycji atomów.
        """
        axis = axis / np.linalg.norm(axis)
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
            - wypełnia brakujące reszty,
            - flipuje jeśli zysk HBonds,
            - zapisuje wynik.
        """
        fixer = PDBFixer(filename=self.input_pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        structure = fixer.topology
        positions = fixer.positions

        for chain in structure.chains():
            for residue in chain.residues():
                resname = residue.name.upper()
                resid = residue.id

                if resname not in ['ASN', 'GLN', 'HIS']:
                    continue

                hbonds_before = self.count_local_hbonds(residue, structure, positions)
                flipped_positions, did_flip = self.flip_and_count_hbonds(resname, residue, structure, positions)

                if not did_flip:
                    continue

                hbonds_after = self.count_local_hbonds(residue, structure, flipped_positions)

                if hbonds_after > hbonds_before:
                    print(f"Zastosowano flip {resname} {resid}, zwiększono lokalne wiązania wodorowe: {hbonds_before} → {hbonds_after}")
                    positions[:] = flipped_positions[:]
                else:
                    print(f"Nie zastosowano flipa {resname} {resid}, brak poprawy: {hbonds_before} → {hbonds_after}")

        output_path = os.path.splitext(self.input_pdb_path)[0] + "_flipped_openmm.pdb"
        with open(output_path, 'w') as out_file:
            PDBFile.writeFile(structure, positions, out_file)
        print(f"Zapisano wynikowy plik do: {output_path}")

# Użycie
# flipper = ResidueFlipper("sciezka/do/pliku.pdb")
# flipper.run()
