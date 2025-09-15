#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import time
from tqdm import tqdm
from shutil import copyfile
from modeller import *
from modeller import environ, model
from modeller.automodel import AutoModel, assess
from modeller.optimizers import molecular_dynamics, actions, conjugate_gradients
from modeller.scripts import complete_pdb
from openbabel import openbabel as ob
import flip_res


# Klasa dziedzicząca po AutoModel – umożliwia późniejszy wybór atomów do modelowania
# Pozwala wskazać, które atomy mają być modelowane
class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self) # domyślnie wybiera wszystkie atomy

# Klasa główna do przetwarzania i naprawy plików PDB
class PDBModelOptimization:
    def __init__(self, project_root):
        """
        Inicjalizacja ścieżek i środowiska Modeller.
        Tworzy foldery robocze, jeśli ich nie ma.
        """

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


    def cleanup_working_files(self, pdb_code):
        """
        Usuwa tymczasowe pliki z katalogu roboczego po zakończeniu pracy.
        """
        for file in os.listdir(self.work_path):
            if file.startswith(pdb_code):
                file_path = os.path.join(self.work_path, file)
                try:
                    os.remove(file_path)
                    print(f"Usunięto pliki robocze: {pdb_code}")
                except Exception as e:
                    print(f"Error przy usuwaniu {file}: {e}")


    def rename_flipped_files(self, pdb_name):
        """
        Zmienia nazwy plików z _flipped_openmm.pdb na .pdb
        """
        flipped_file = os.path.join(self.output_path, f"{pdb_name}_flipped_openmm.pdb")
        target_file = os.path.join(self.output_path, f"{pdb_name}.pdb")

        if os.path.exists(flipped_file):
            if os.path.exists(target_file):
                os.remove(target_file)
            os.rename(flipped_file, target_file)
            print(f"Zmieniono nazwę: {flipped_file} -> {target_file}")


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
        self.aln = aln


    def add_hydrogens(self, pdb_path):
        """
        Dodaje wodory do pliku PDB, próbując najpierw przez subprocess (obabel CLI),
        a potem przez Open Babel Python API. Nadpisuje plik wejściowy.
        """
        try:
            subprocess.run(["obabel", pdb_path, "-O", pdb_path, "-h"], check=True)
            print(f"Dodano wodory za pomocą Open Babel CLI: {pdb_path}")
            return
        except subprocess.CalledProcessError as e:
            print("Błąd Open Babel CLI:", e)
        except Exception as e:
            print("Nieoczekiwany błąd podczas uruchamiania obabel:", e)

        # Fallback – OpenBabel Python API
        try:
            obConversion = ob.OBConversion()
            obConversion.SetInAndOutFormats("pdb", "pdb")

            mol = ob.OBMol()
            if not obConversion.ReadFile(mol, pdb_path):
                print(f"Błąd przy wczytywaniu pliku: {pdb_path}")
                return

            mol.AddHydrogens()

            if obConversion.WriteFile(mol, pdb_path):
                print(f"Dodano wodory za pomocą OpenBabel Python API: {pdb_path}")
            else:
                print(f"Błąd przy zapisie pliku: {pdb_path}")
        except Exception as e:
            print(f"Błąd OpenBabel Python API: {e}")


    def remove_duplicate_atoms(self, pdb_path):
        """
        Filtruje plik PDB z powielonych atomów.
        Unika błędów przy dalszym przetwarzaniu.
        """
        seen = set()
        cleaned_lines = []

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_name = line[12:16].strip()
                    res_id = line[17:26]  # Łańcuch + numer + insertion code
                    key = (res_id, atom_name)

                    if key in seen:
                        continue  # pomiń duplikat
                    seen.add(key)

                cleaned_lines.append(line)

        with open(pdb_path, 'w') as f:
            f.writelines(cleaned_lines)

        print(f"Usunięto duplikaty atomów w: {pdb_path}")

    def optimize_heavy_atom(self, model, code):
        """
        Optymalizuje geometrię atomów ciężkich (bez wodorów) za pomocą gradientów sprzężonych.
        """
        heavy_atoms = selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])
        model.restraints.make(heavy_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        cg.optimize(heavy_atoms, max_iterations=30)
        final_path = os.path.join(self.output_path, code + ".pdb")
        model.write(file=final_path)
        print(f"Zoptymalizowano atomy ciężkie i zapisano: {final_path}")


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
        cg.optimize(all_atoms, max_iterations=300, actions=[actions.trace(10, trace_file)])

        # Dynamika molekularna z rampowaniem temperatury
        for temp in [100, 200, 300]:
            md.optimize(all_atoms, temperature=temp, max_iterations=200,
                        actions=[actions.write_structure(20, code + ".dyn%04d.pdb"),
                                 actions.trace(20, trace_file)])

        # Ponowna optymalizacja
        cg.optimize(all_atoms, max_iterations=200, actions=[actions.trace(10, trace_file)])


        final_path = os.path.join(self.output_path, code + ".pdb")
        model.write(file=final_path)
        print(f"Zoptymalizowano pełną strukturę i zapisano: {final_path}")


    def fill_missing_residues_and_atoms(self, single_file=None):
        """
        Główna funkcja przetwarzająca wszystkie pliki PDB:
             - uzupełnia brakujące reszty,
             - dodaje brakujące atomy,
             - dodaje atomy wodoru,
             - flipuje reszty HIS/ASN/GLN,
             - przeprowadza optymalizację struktury.
        """
        files = [single_file] if single_file else os.listdir(self.input_path)
        stats = {}

        for filename in tqdm(files):
            if not filename.endswith(".pdb"):
                continue

            print(f"Przetwarzanie pliku: {filename}")
            start = time.perf_counter()

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

            prepare_alignment_end = time.perf_counter()
            prepare_alignment_time = prepare_alignment_end - start
            print("prepare_alignment_time trwał: ", prepare_alignment_time)

            # Zapisanie pliku wyjsciowego
            model_output_file = f"{temple_name}.B99990001.pdb"
            if os.path.isfile(model_output_file):
                final_output_file = os.path.join(self.output_path, f"{pdb_name}.pdb")
                if os.path.exists(final_output_file):
                    os.remove(final_output_file)
                os.rename(model_output_file, final_output_file)
                print(f"Uzupełniony model zapisany do: {final_output_file}")

                # 1) Dodanie wodorow (tylko jeśli ich nie ma)
                with open(final_output_file, 'r') as f:
                    lines = f.readlines()

                has_hydrogens = any(
                    line.startswith("ATOM") and line[76:78].strip() == "H"
                    for line in lines
                )

                # Sprawdzenie, czy są wodory — jeśli nie, to je dodaj
                if not has_hydrogens:
                    self.add_hydrogens(final_output_file)
                    self.remove_duplicate_atoms(final_output_file)
                else:
                    print(f"{final_output_file} już zawiera atomy wodoru – pomijam dodawanie.")
                    self.remove_duplicate_atoms(final_output_file)

                code = os.path.splitext(os.path.basename(final_output_file))[0]
                self.env.io.atom_files_directory.append(os.path.dirname(final_output_file))

                hydrogens_added_end = time.perf_counter()
                hydrogens_added_time = hydrogens_added_end - prepare_alignment_end
                print("hydrogens_added_time trwał: ", hydrogens_added_time)

                # 2) Wczytanie uzupelnionej struktury
                model = complete_pdb(self.env, code)

                flip_res.ResidueFlipper(final_output_file, self.aln)

                model_loaded_end = time.perf_counter()
                model_loaded_time = model_loaded_end - hydrogens_added_end
                print("model_loaded_time trwał: ", model_loaded_time)

                # 3) Flipy odpowiednich reszt aminokwasowych
                flipper = flip_res.ResidueFlipper(final_output_file, self.aln)
                flipper.run()

                residues_flipped_end = time.perf_counter()
                residues_flipped_time = residues_flipped_end - model_loaded_end
                print("residues_flipped_time trwał: ", residues_flipped_time)

                # 4) Optymalizacje
                self.optimize_heavy_atom(model, code)
                self.optimize_full_structure(model, code)

                optimized_end = time.perf_counter()
                optimized_time = optimized_end - residues_flipped_end

                # 5) Zmiana nazwy plików na poprawne (podmienia pliki)
                self.rename_flipped_files(pdb_name)

                # 6) Czyszczenie plików roboczych
                self.cleanup_working_files(pdb_name)

                # 7) Licznik
                files_cleaned_end = time.perf_counter()
                files_cleaned_time = files_cleaned_end - optimized_end
                print("files_cleaned_time trwał: ", files_cleaned_time)

                file_stats = {
                    "prepare_alignment_time": prepare_alignment_time,
                    "hydrogens_added_time": hydrogens_added_time,
                    "model_loaded_time": model_loaded_time,
                    "residues_flipped_time": residues_flipped_time,
                    "optimized_time": optimized_time,
                    "files_cleaned_time": files_cleaned_time
                }
                print(file_stats)
                stats[filename] = file_stats
                end = time.perf_counter()
                print("całość: ", end - start)
            else:
                print(f"Błąd: nie znaleziono modelu dla {pdb_name}!")

        print(stats)



if __name__ == '__main__':
    # Parser argumentów: pozwala uruchomić skrypt dla 1 pliku lub wszystkich.
    parser = argparse.ArgumentParser(description="Naprawianie struktury PDB (reszty, wodory, flipy)")
    parser.add_argument("-f", "--file", type=str, help="Nazwa pliku .pdb do przetworzenia")
    args = parser.parse_args()

    project_root = os.getcwd()
    processor = PDBModelOptimization(project_root)

    if args.file:
        processor.fill_missing_residues_and_atoms(single_file=args.file)
    else:
        processor.fill_missing_residues_and_atoms()