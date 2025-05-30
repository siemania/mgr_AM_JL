import os
from shutil import copyfile
from modeller import *
from modeller.automodel import AutoModel, assess


class MyModel(AutoModel):
    def select_atoms(self):
        return selection(self)


def prepare_alignment(env, pdb_code, temple_code):
    """Tworzy alignment zawierający model źródłowy i docelowy"""
    mdl = model(env, file=pdb_code)
    aln = alignment(env)

    # Dodanie znanej struktury do alignment
    aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code)

    # Dodanie  struktury do poprawienia do alignment
    aln.append_model(mdl, atom_files=pdb_code, align_codes=temple_code)

    # Zapis alignment w formacie .ali
    aln.write(file='alignment.ali')
    print("alignment.ali zostal zapisany")


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
        temple_name = f"{pdb_name}_fill" #plik szablonu

        env.io.atom_files_directory = ['.', os.path.dirname(input_path)]

        os.chdir(work_path)

        if not os.path.exists(filename):
            copyfile(pdb_path, filename)

        # Przygotowanie alignment
        prepare_alignment(env, pdb_name, temple_name)

        # Modelowanie
        a = MyModel(env,
                    alnfile='alignment.ali',
                    knowns=pdb_name,
                    sequence=temple_name,
                    assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        #Zapisanie pliku wyjsciowego
        model_output_file = f"{temple_name}.B99990001.pdb"
        if os.path.isfile(model_output_file):
            final_output_file = os.path.join(output_path, f"{pdb_name}.pdb")
            if os.path.exists(final_output_file):
                os.remove(final_output_file)
            os.rename(model_output_file, final_output_file)
            print(f"Uzupełniony model zapisany do: {final_output_file}")
        else:
            print(f"Błąd: nie znaleziono modelu dla {pdb_name}!")
if __name__ == '__main__':
    # Ścieżki plików
    username = os.environ.get('USERNAME') or os.environ.get('USER')
    PROJECT_ROOT = rf"C:\Users\{username}\PycharmProjects\mgr_AM_JL"
    input_pdb = os.path.join(PROJECT_ROOT, 'pdb_files')
    output_pdb = os.path.join(PROJECT_ROOT, 'fixed_pdb')


    fill_missing_residues_and_atoms(input_pdb, output_pdb, PROJECT_ROOT)
