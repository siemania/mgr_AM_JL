import os
from shutil import copyfile
from modeller import *
from modeller import environ, model
from modeller.automodel import AutoModel, assess
from modeller.scripts import complete_pdb
from modeller.optimizers import molecular_dynamics, actions, conjugate_gradients
from Docking import flip_HIS



class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self)

class PDBModelOptimization:
    def __init__(self, project_root):
        self.project_root = project_root
        self.input_path = os.path.join(project_root, 'pdb_files') # Folder z wejsciowymi plikami PDB
        self.output_path = os.path.join(project_root, 'fixed_pdb') # Folder na naprawione pliki PDB
        self.work_path = os.path.join(project_root, 'work_folder') # Folder roboczy

        os.makedirs(self.input_path, exist_ok=True)
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.work_path, exist_ok=True)

        # Inicjalizacja środowiska Modeller
        log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
        self.env = environ()
        self.env.io.atom_files_directory = ['.', self.input_path]
        self.env.libs.topology.read(file='$(LIB)/top_heav.lib')
        self.env.libs.parameters.read(file='$(LIB)/par.lib')
        self.env.edat.dynamic_sphere = True

    def prepare_alignment(self, env, pdb_code, temple_code):
        mdl = model(env, file=pdb_code)
        aln = alignment(env)
        aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code)
        aln.append_model(mdl, atom_files=pdb_code, align_codes=temple_code)
        aln.write(file='alignment.ali')
        print("alignment.ali zostal zapisany")

    def swap_atom_coords(self, atom1, atom2):
        """Zamienia wspolrzedne dwoch podanych atomow."""
        atom1.x, atom2.x = atom2.x, atom1.x
        atom1.y, atom2.y = atom2.y, atom1.y
        atom1.z, atom2.z = atom2.z, atom1.z

    def flip_residues (self, residue):
        """Flipuje reszty HIS, ASN i GLN przez zamianę współrzędnych wybranych atomów."""
        atoms = residue.atom_dict

        if residue.name == 'ASN':
            # Zamienia wspolrzedne OD1 i OD2
            if 'OD1' in atoms and 'OD2' in atoms:
                self.swap_atom_coords(self, atoms['OD1'], atoms['OD2'])
                print(f"Wykonano flip ASN: {residue.num}")

        elif residue.name == 'GLN':
            # Zamienia wspolrzedne OE1 i OE2
            if 'OE1' in atoms and 'OE2' in atoms:
                self.swap_atom_coords(self, atoms['OE1'], atoms['OE2'])
                print(f"Wykonano flip GLN: {residue.num}")

    def optimize_heavy_atom(self, model, code):
        heavy_atoms = selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])
        model.restraints.make(heavy_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        cg.optimize(heavy_atoms, max_iterations=30)
        model.write(file=code + '_optHeavy.pdb')
        print(f"Zoptymalizowano atomy ciężkie: {code}_optHeavy.pdb")


    def optimize_full_structure(self, model, code):

        model.write(file=code + ".ini")

        all_atoms = Selection(model)
        model.restraints.make(all_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(code + ".full_opt.log", "w")

        cg.optimize(all_atoms, max_iterations=30, actions=[actions.trace(5, trace_file)])
        md.optimize(all_atoms, temperature=300, max_iterations=50,
                    actions=[actions.write_structure(10, code + ".dyn%04d.pdb"),
                             actions.trace(10, trace_file)])
        cg.optimize(all_atoms, max_iterations=20, actions=[actions.trace(5, trace_file)])

        model.write(file=code + '_opt_full.pdb')
        print(f"Zoptymalizowano całą strukturę: {code}_opt_full.pdb")


    def fill_missing_residues_and_atoms(self):

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


                code = pdb_name
                model = complete_pdb(self.env, code)
                # Flipy aminokwasow
                for res in model.residues:
                    if residue.name in ['ASN', 'GLN']:
                        self.flip_residues(res)
                    elif residue.name == 'HIS':
                        flip_HIS(self.input_path, self.work_path, res)
                # Optymalizacje
                self.optimize_heavy_atom(model, code)
                self.optimize_full_structure(model, code)
            else:
                print(f"Błąd: nie znaleziono modelu dla {pdb_name}!")


if __name__ == '__main__':

    project_root = os.getcwd()
    processor = PDBModelOptimization(project_root)
    processor.fill_missing_residues_and_atoms()
