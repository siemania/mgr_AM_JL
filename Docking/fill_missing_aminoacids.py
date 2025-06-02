import os
from shutil import copyfile
from modeller import *
from modeller.automodel import AutoModel, assess
from modeller.scripts import complete_pdb
from modeller.optimizers import conjugate_gradients, molecular_dynamics, actions, ConjugateGradients


class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self)


def prepare_alignment(env, pdb_code, temple_code):
    mdl = model(env, file=pdb_code)
    aln = alignment(env)
    aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code)
    aln.append_model(mdl, atom_files=pdb_code, align_codes=temple_code)
    aln.write(file='alignment.ali')
    print("alignment.ali zostal zapisany")


def optimize_with_heavy_atom_restraints(pdb_file_path):
    env = environ()
    env.io.atom_files_directory = ['.']
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    env.edat.dynamic_sphere = True

    code = os.path.splitext(os.path.basename(pdb_file_path))[0]
    mdl = complete_pdb(env, code)

    heavy_atoms = selection(*[atom for atom in mdl.atoms if not atom.name.strip().startswith('H')])
    mdl.restraints.make(heavy_atoms, restraint_type='stereo', spline_on_site=False)

    cg = conjugate_gradients(output='REPORT')
    cg.optimize(heavy_atoms, max_iterations=30)
    mdl.write(file=code + '_optHeavy.pdb')
    print(f"Zoptymalizowano atomy ciężkie: {code}_optHeavy.pdb")


def optimize_full_structure(pdb_file_path):
    log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
    env = environ()
    env.io.atom_files_directory = ['.']
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    code = os.path.splitext(os.path.basename(pdb_file_path))[0]
    mdl = complete_pdb(env, code)
    mdl.write(file=code + ".ini")

    all_atoms = Selection(mdl)
    mdl.restraints.make(all_atoms, restraint_type='stereo', spline_on_site=False)

    cg = conjugate_gradients(output='REPORT')
    md = molecular_dynamics(output='REPORT')
    trace_file = open(code + ".full_opt.log", "w")

    cg.optimize(all_atoms, max_iterations=30, actions=[actions.trace(5, trace_file)])
    md.optimize(all_atoms, temperature=300, max_iterations=50,
                actions=[actions.write_structure(10, code + ".dyn%04d.pdb"),
                         actions.trace(10, trace_file)])
    cg.optimize(all_atoms, max_iterations=20, actions=[actions.trace(5, trace_file)])

    mdl.write(file=code + '_opt_full.pdb')
    print(f"Zoptymalizowano całą strukturę: {code}_opt_full.pdb")


def fill_missing_residues_and_atoms(input_path, output_path, project_root):
    log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
    env = environ()

    work_path = os.path.join(project_root, 'work_folder')
    os.makedirs(work_path, exist_ok=True)
    os.makedirs(output_path, exist_ok=True)

    for filename in os.listdir(input_path):
        if not filename.endswith(".pdb"):
            continue

        print(f"Przetwarzanie pliku: {filename}")

        pdb_path = os.path.join(input_path, filename)
        pdb_name = os.path.splitext(filename)[0]
        temple_name = f"{pdb_name}_fill"

        env.io.atom_files_directory = ['.', os.path.dirname(input_path)]

        os.chdir(work_path)

        if not os.path.exists(filename):
            copyfile(pdb_path, filename)

        prepare_alignment(env, pdb_name, temple_name)

        a = MyModel(env,
                    alnfile='alignment.ali',
                    knowns=pdb_name,
                    sequence=temple_name,
                    assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        model_output_file = f"{temple_name}.B99990001.pdb"
        if os.path.isfile(model_output_file):
            final_output_file = os.path.join(output_path, f"{pdb_name}.pdb")
            if os.path.exists(final_output_file):
                os.remove(final_output_file)
            os.rename(model_output_file, final_output_file)
            print(f"Uzupełniony model zapisany do: {final_output_file}")

            #Optymalizacje
            optimize_with_heavy_atom_restraints(final_output_file)
            optimize_full_structure(final_output_file)
        else:
            print(f"Błąd: nie znaleziono modelu dla {pdb_name}!")


if __name__ == '__main__':
    username = os.environ.get('USERNAME') or os.environ.get('USER')
    PROJECT_ROOT = rf"C:\\Users\\{username}\\PycharmProjects\\mgr_AM_JL"
    input_pdb = os.path.join(PROJECT_ROOT, 'pdb_files')
    output_pdb = os.path.join(PROJECT_ROOT, 'fixed_pdb')

    fill_missing_residues_and_atoms(input_pdb, output_pdb, PROJECT_ROOT)
