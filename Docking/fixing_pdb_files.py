import os
from shutil import copyfile
import numpy as np
from modeller import *
from modeller import environ, model
from modeller.automodel import AutoModel, assess
from modeller.scripts import complete_pdb
from modeller.optimizers import molecular_dynamics, actions, conjugate_gradients
from openbabel import openbabel
from openmm.app import PDBFile
from openmm import unit
from pdbfixer import PDBFixer



# Klasa dziedzicząca po AutoModel – umożliwia późniejszy wybór atomów do modelowania
class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self) # domyślnie wybiera wszystkie atomy

class PDBModelOptimization:
    def __init__(self, project_root):
        # Ścieżki do folderów wejściowych, wyjściowych i roboczych
        self.project_root = project_root
        self.input_path = os.path.join(project_root, 'pdb_files') # Folder z wejsciowymi plikami PDB
        self.output_path = os.path.join(project_root, 'fixed_pdb') # Folder na naprawione pliki PDB
        self.work_path = os.path.join(project_root, 'work_folder') # Folder roboczy

        os.makedirs(self.input_path, exist_ok=True)
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.work_path, exist_ok=True)

        # Inicjalizacja środowiska Modeller
        log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
        self.env = environ()
        self.env.io.atom_files_directory = ['.', self.input_path]
        self.env.libs.topology.read(file='$(LIB)/top_heav.lib')
        self.env.libs.parameters.read(file='$(LIB)/par.lib')
        self.env.edat.dynamic_sphere = True

    def prepare_alignment(self, env, pdb_code, temple_code):
        """
        Tworzy alignment pomiędzy znaną strukturą (template) a sekwencją docelową (target).
        Potrzebne do działania AutoModel.
        """
        mdl = model(env, file=pdb_code)
        aln = alignment(env)
        aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code)
        aln.append_model(mdl, atom_files=pdb_code, align_codes=temple_code)
        aln.write(file='alignment.ali')
        print("alignment.ali zostal zapisany")


    def add_hydrogens(self, pdb_path):
        output_path = pdb_path.replace(".pdb", "_withH.pdb")

        try:
            fixer = PDBFixer(filename=pdb_path)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH=7.0)

            with open(output_path, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

            print(f"Dodano wodory za pomocą PDBFixer: {output_path}")
            return output_path

        except Exception as e:
            print(f"Błąd przy dodawaniu wodorów przez PDBFixer: {e}")
            return pdb_path

    def rotate_atoms(self, positions, center, axis, angle_rad):
        """Obraca pozycje wokół zadanej osi i punktu."""
        # Zamień center i axis na ndarray
        axis = np.array(axis)
        center = np.array(center)

        # Normalizacja osi
        axis = axis / np.linalg.norm(axis)

        rotated = []

        for pos in positions:
            rel = pos - center
            rot = (
                    np.cos(angle_rad) * rel +
                    np.sin(angle_rad) * np.cross(axis, rel) +
                    (1 - np.cos(angle_rad)) * np.dot(axis, rel) * axis
            )
            rotated.append(rot + center)

        return [unit.Quantity(r, unit.nanometer) for r in rotated]

    def flip_HIS(self, input_file):
        try:
            fixer = PDBFixer(filename=input_file)
            structure = fixer.topology
            positions = fixer.positions

            changed = False
            for residue in structure.residues():
                if residue.name not in ('HIS', 'HID', 'HIE', 'HIP'):
                    continue

                atom_dict = {atom.name: atom for atom in residue.atoms()}
                required_atoms = ['CG', 'CD2', 'ND1', 'CE1', 'NE2']

                if not all(name in atom_dict for name in required_atoms):
                    continue  # pomiń niekompletne reszty

                cg_idx = list(structure.atoms()).index(atom_dict['CG'])
                cd2_idx = list(structure.atoms()).index(atom_dict['CD2'])

                axis = positions[cd2_idx].value_in_unit(unit.nanometer) - positions[cg_idx].value_in_unit(
                    unit.nanometer)
                center = positions[cg_idx].value_in_unit(unit.nanometer)

                ring_atom_indices = [list(structure.atoms()).index(atom_dict[name]) for name in
                                     ['ND1', 'CD2', 'CE1', 'NE2']]
                ring_positions = [positions[i].value_in_unit(unit.nanometer) for i in ring_atom_indices]


                new_positions = self.rotate_atoms(ring_positions, center, axis, np.pi)

                for i, idx in enumerate(ring_atom_indices):
                    positions[idx] = new_positions[i]

                changed = True

            if changed:
                output_file = input_file.replace(".pdb", "_rotatedHIS.pdb")
                with open(output_file, 'w') as f:
                    PDBFile.writeFile(structure, positions, f)
                return output_file
            else:
                return input_file #brak zmian, zwraca orginał
        except Exception as e:
            print(f"Błąd przy obracaniu HIS: {e}")
            return input_file

    def swap_atom_coords(self, atom1, atom2):
        """
        Zamienia wspolrzedne (x, y, z) dwoch podanych atomow.
        """
        atom1.x, atom2.x = atom2.x, atom1.x
        atom1.y, atom2.y = atom2.y, atom1.y
        atom1.z, atom2.z = atom2.z, atom1.z

    def flip_residues (self, residue, pdb_path=None):
        """
        Flipuje reszty GLN, ASN i HIS przez zamianę współrzędnych wybranych atomów.
        Pomaga ustawić poprawną orientację grup funkcyjnych.
        Zachodzą zmiany miedzy atomami:
            - ASN: OD1 <-> ND2
            - GLN: OE1 <-> NE2
            - HIS:  ND1 <-> NE2 (uproszczony sposob)
         """
        atom_dict = {a.name.strip(): a for a in residue.atoms}
        name = residue.name.strip().upper()

        if name == 'ASN':

            if 'OD1' in atom_dict and 'ND2' in atom_dict:
                self.swap_atom_coords(atom_dict['OD1'], atom_dict['ND2'])
                print(f"Wykonano flip ASN: {residue.num}")

        elif name == 'GLN':

            if 'OE1' in atom_dict and 'NE2' in atom_dict:
                self.swap_atom_coords(atom_dict['OE1'], atom_dict['NE2'])
                print(f"Wykonano flip GLN: {residue.num}")
        elif name == 'HIS' and pdb_path:
            print(f"Obracanie HIS w {residue.num} przez PDBFixer")
            rotated_path = self.flip_HIS(pdb_path)
            if rotated_path:
                print(f"Obrócono pierścień histydyny: {residue.num}")

    def optimize_heavy_atom(self, model, code):
        """
        Optymalizuje geometrię atomów ciężkich (bez wodorów) za pomocą gradientów sprzężonych.
        """
        heavy_atoms = selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])
        model.restraints.make(heavy_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        cg.optimize(heavy_atoms, max_iterations=30)
        model.write(file=code + '_optHeavy.pdb')
        print(f"Zoptymalizowano atomy ciężkie: {code}_optHeavy.pdb")


    def optimize_full_structure(self, model, code):
        """
        Pełna optymalizacja struktury: gradienty sprzężone + dynamiczne optymalizacje molekularne.
        """
        model.write(file=code + ".ini")

        all_atoms = Selection(model)
        model.restraints.make(all_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(code + ".full_opt.log", "w")

        # Gradienty sprzężone
        cg.optimize(all_atoms, max_iterations=30, actions=[actions.trace(5, trace_file)])

        # Dynamika molekularna
        md.optimize(all_atoms, temperature=300, max_iterations=50,
                    actions=[actions.write_structure(10, code + ".dyn%04d.pdb"),
                             actions.trace(10, trace_file)])

        # Ponowna optymalizacja
        cg.optimize(all_atoms, max_iterations=20, actions=[actions.trace(5, trace_file)])

        model.write(file=code + '_opt_full.pdb')
        print(f"Zoptymalizowano całą strukturę: {code}_opt_full.pdb")


    def fill_missing_residues_and_atoms(self):
        """
        Główna funkcja przetwarzająca wszystkie pliki PDB:
             - uzupełnia brakujące reszty,
             - dodaje brakujące atomy,
             - dodaje atomy wodoru,
             - flipuje reszty HIS/ASN/GLN,
             - przeprowadza optymalizację struktury.
        """
        for filename in os.listdir(self.input_path):
            if not filename.endswith(".pdb"):
                continue

            print(f"Przetwarzanie pliku: {filename}")

            pdb_path = os.path.join(self.input_path, filename)
            pdb_name = os.path.splitext(filename)[0]
            temple_name = f"{pdb_name}_fill"

            self.env.io.atom_files_directory = ['.', self.input_path]

            os.chdir(self.work_path)

            if not os.path.exists(filename):
                copyfile(pdb_path, filename)

            self.prepare_alignment(self.env, pdb_name, temple_name)

            a = MyModel(self.env,
                        alnfile='alignment.ali',
                        knowns=pdb_name,
                        sequence=temple_name,
                        assess_methods=(assess.DOPE, assess.GA341))
            a.starting_model = 1
            a.ending_model = 1
            a.make()

            # Zapisanie pliku wyjsciowego
            model_output_file = f"{temple_name}.B99990001.pdb"
            if os.path.isfile(model_output_file):
                final_output_file = os.path.join(self.output_path, f"{pdb_name}.pdb")
                if os.path.exists(final_output_file):
                    os.remove(final_output_file)
                os.rename(model_output_file, final_output_file)
                print(f"Uzupełniony model zapisany do: {final_output_file}")


                # 1) Dodanie wodorow
                pdb_path_with_H = self.add_hydrogens(final_output_file)

                code = os.path.splitext(os.path.basename(pdb_path_with_H))[0]
                self.env.io.atom_files_directory.append(os.path.dirname(pdb_path_with_H))

                # 2) Wczytanie uzupelnionej struktury
                model = complete_pdb(self.env, code)

                # 3) Flipy odpowiednich reszt aminokwasowych
                for res in model.residues:
                    self.flip_residues(res, pdb_path_with_H)

                # 4) Optymalizacje
                self.optimize_heavy_atom(model, code)
                self.optimize_full_structure(model, code)
            else:
                print(f"Błąd: nie znaleziono modelu dla {pdb_name}!")


if __name__ == '__main__':

    project_root = os.getcwd()
    processor = PDBModelOptimization(project_root)
    processor.fill_missing_residues_and_atoms()
