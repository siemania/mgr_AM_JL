#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import time
from tqdm import tqdm
from shutil import copyfile
from modeller import *
from modeller import environ, model
from modeller.automodel import AutoModel, assess
from modeller.optimizers import molecular_dynamics, actions, conjugate_gradients
from modeller.scripts import complete_pdb
from scipy.spatial import cKDTree
import numpy as np
import flip_res


# Klasa MyModel dziedziczy po AutoModel (Modeller)
# Umożliwia selekcję konkretnych atomów do optymalizacji (domyślnie wybiera wszystkie)
class MyModel(AutoModel):
    def select_atoms(self):
        return Selection(self) # domyślnie wybiera wszystkie atomy

# Klasa główna do przetwarzania i naprawy plików PDB
class PDBModelOptimization:
    """
        Klasa do automatycznego przetwarzania, uzupełniania i optymalizacji struktur PDB.
        Proces obejmuje:
          - uzupełnianie brakujących reszt i atomów (Modeller)
          - dodawanie wodorów (top_allh.lib)
          - flipowanie reszt HIS/ASN/GLN (analiza orientacji grup amidowych)
          - optymalizacje geometryczne i dynamiczne (conjugate gradients, MD)
          - naprawę i zapis finalnych struktur PDB
    """
    def __init__(self, project_root):
        """
            Inicjalizacja środowiska Modeller oraz konfiguracja ścieżek roboczych.
            Tworzy wymagane foldery robocze, jeśli nie istnieją.
        """

        # Ścieżki do folderów wejściowych, wyjściowych i roboczych
        self.project_root = project_root
        self.input_path = os.path.join(project_root, 'pdb_files') # Folder z wejsciowymi plikami PDB
        self.output_path = os.path.join(project_root, 'pdb_files') # Folder na naprawione pliki PDB
        self.work_path = os.path.join(project_root, 'work_folder') # Folder roboczy

        os.makedirs(self.input_path, exist_ok=True)
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.work_path, exist_ok=True)

        # Inicjalizacja środowiska Modeller
        log.level(output=0, notes=0, warnings=0, errors=1, memory=0)
        self.env = environ()
        self.env.io.atom_files_directory = ['.', self.input_path]
        self.env.libs.topology.read(file='$(LIB)/top_allh.lib')
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
            Zastępuje nazwę pliku *_flipped_openmm.pdb standardową nazwą *.pdb.
            Pozwala utrzymać spójność nazw po flipowaniu reszt aminokwasowych.
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
            Tworzy alignment pomiędzy strukturą wzorcową (template) a sekwencją docelową (target).
            Alignment zapisywany jest w pliku 'alignment.ali' i wykorzystywany przez AutoModel.
        """
        mdl = model(env, file=pdb_code)
        aln = alignment(env)

        aln.append_model(mdl, atom_files=pdb_code, align_codes=pdb_code)
        aln.append_model(mdl, atom_files=pdb_code, align_codes=temple_code)

        aln.write(file='alignment.ali')

        print("alignment.ali zostal zapisany")
        self.aln = aln





    def optimize_heavy_atom(self, model, code):
        """
            Przeprowadza lokalną optymalizację geometrii atomów ciężkich (bez H)
            przy użyciu metody gradientów sprzężonych (Conjugate Gradients).
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
            Wykonuje pełną optymalizację modelu:
                - relaksacja struktury za pomocą gradientów sprzężonych,
                - symulacja dynamiki molekularnej (MD) z rampowaniem temperatury,
                - ponowna optymalizacja po MD.
            Proces minimalizuje energię układu, poprawiając geometrię lokalną i globalną.
        """
        model.write(file=code + ".ini")

        all_atoms = Selection(model)
        model.restraints.make(all_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(code + ".full_opt.log", "w")

        # Gradienty sprzężone
        cg.optimize(all_atoms, max_iterations=50, actions=[actions.trace(10, trace_file)])

        # Dynamika molekularna z rampowaniem temperatury
        for temp in [150, 300]:
            md.optimize(all_atoms, temperature=temp, max_iterations=100,
                        actions=[actions.write_structure(20, code + ".dyn%04d.pdb"),
                                 actions.trace(20, trace_file)])

        # Ponowna optymalizacja
        cg.optimize(all_atoms, max_iterations=50, actions=[actions.trace(10, trace_file)])


        final_path = os.path.join(self.output_path, code + ".pdb")
        model.write(file=final_path)
        print(f"Zoptymalizowano pełną strukturę i zapisano: {final_path}")

    def find_active_site(self, model, ligand_pdb_path, cutoff=5.0):
        """
            Identyfikuje reszty centrum aktywnego na podstawie odległości od liganda.
            Wykorzystuje drzewo KDTree (scipy.spatial) do efektywnego wyszukiwania sąsiadujących reszt.

            Parametry:
                model : obiekt Modeller model
                    Model białka
                ligand_pdb_path : str
                    Ścieżka do pliku liganda (PDB)
                cutoff : float
                    Promień (Å), w którym wyszukiwane są reszty w centrum aktywnym
        """

        ligand_atoms = [a for r in model.residues for a in r.atoms]

        if not ligand_atoms:
            print("Nie znaleziono atomów liganda w pliku:", ligand_pdb_path)
            return Selection()

        # współrzędne atomów liganda
        ligand_coords = np.array([[a.x, a.y, a.z] for a in ligand_atoms])

        # współrzędne atomów modelu (bez wody)
        model_atoms = [a for a in model.atoms if a.residue.name != "HOH"]
        model_coords = np.array([[a.x, a.y, a.z] for a in model_atoms])

        # budowa KDTree dla modelu
        tree = cKDTree(model_coords)


        nearby_indices = tree.query_ball_point(ligand_coords, cutoff)
        nearby_residues = {model_atoms[i].residue for idx_list in nearby_indices for i in idx_list}

        active_site = Selection(*nearby_residues)
        print(f"Znaleziono {len(nearby_residues)} reszt w centrum aktywnym (KDTree optymalizacja)")
        return active_site


    def optimize_with_restraints(self, model, code, pdb_name):
        """
            Optymalizacja struktury z nałożonymi więzami:
                - szkielet białka (backbone) ograniczony,
                - centrum aktywne zamrożone (atomy liganda i najbliższe reszty),
                - pozostałe elementy poddane relaksacji.

            Celem jest poprawa lokalnej geometrii bez zaburzania konformacji funkcjonalnych.
        """

        all_atoms = Selection(model)

        # Backbone (N, CA, C, O)
        backbone = Selection(
            *[atom for atom in model.atoms if atom.name.strip() in ("N", "CA", "C", "O")]
        )

        # jeśli masz ligand w folderze /ligands
        ligand_path = os.path.join(self.project_root, "ligands", f"{pdb_name}_ligand.pdb")
        active_site = self.find_active_site(model, ligand_path)


        # Więzy: backbone + stereo restraints
        model.restraints.make(backbone, restraint_type='stereo', spline_on_site=False)


        # Zamrożenie centrum aktywnego
        for atom in active_site:
            atom.fix = True

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(code + ".restrained_opt.log", "w")

        # Gradienty + MD + Gradienty
        cg.optimize(all_atoms, max_iterations=100, actions=[actions.trace(10, trace_file)])
        md.optimize(all_atoms, temperature=300, max_iterations=300,
                    actions=[actions.write_structure(50, code + ".dyn%04d.pdb"),
                             actions.trace(20, trace_file)])
        cg.optimize(all_atoms, max_iterations=200, actions=[actions.trace(10, trace_file)])

        final_path = os.path.join(self.output_path, code + "_restrained.pdb")
        model.write(file=final_path)
        print(f"Zapisano zoptymalizowaną strukturę (z więzami): {final_path}")

    def fill_missing_residues_and_atoms(self, single_file=None):
        """
            Główna funkcja przetwarzająca pojedyncze lub wszystkie pliki PDB.
            Wykonuje kompletny pipeline naprawczy:

                1. Uzupełnienie brakujących reszt/atomów (Modeller)
                2. Dodanie atomów wodoru (top_allh.lib)
                3. Flipowanie reszt ASN, GLN, HIS (flip_res.py)
                4. Optymalizację struktury z więzami (restrained optimization)
                5. Pełną optymalizację energetyczną (full structure optimization)
                6. Sprzątanie plików tymczasowych

            Zwraca statystyki czasów poszczególnych etapów.
        """

        files = [single_file] if single_file else os.listdir(self.input_path)
        stats = {}
        print(self.input_path)
        print(self.output_path)

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

            # 1) Zapisanie pliku wyjsciowego

            model_output_file = f"{temple_name}.B99990001.pdb"
            if os.path.isfile(model_output_file):
                final_output_file = os.path.join(self.output_path, f"{pdb_name}.pdb")
                if os.path.exists(final_output_file):
                    os.remove(final_output_file)
                os.rename(model_output_file, final_output_file)
                print(f"Uzupełniony model zapisany do: {final_output_file}")


                code = os.path.splitext(os.path.basename(final_output_file))[0]
                self.env.io.atom_files_directory.append(os.path.dirname(final_output_file))

                hydrogens_added_end = time.perf_counter()
                hydrogens_added_time = hydrogens_added_end - prepare_alignment_end
                print("hydrogens_added_time trwał: ", hydrogens_added_time)

                # 2) Wczytanie uzupelnionej struktury
                model = complete_pdb(self.env, code)
                model.write(file=f"{self.output_path}/{pdb_name}_withH.pdb")
                print("Modeller automatycznie dodał wodory (top_allh.lib).")



                model_loaded_end = time.perf_counter()
                model_loaded_time = model_loaded_end - hydrogens_added_end
                print("model_loaded_time trwał: ", model_loaded_time)

                # 3) Flipy odpowiednich reszt aminokwasowych
                flipper = flip_res.ResidueFlipper(final_output_file, self.aln)
                flipper.run()

                # Wczytanie ponownie strukturę po flipowaniu

                flipped_file = f"{self.output_path}/{pdb_name}_flipped_openmm.pdb"
                if os.path.exists(flipped_file):
                    model = complete_pdb(self.env, f"{pdb_name}_flipped_openmm")
                else:
                    model = complete_pdb(self.env, f"{pdb_name}_withH")

                residues_flipped_end = time.perf_counter()
                residues_flipped_time = residues_flipped_end - model_loaded_end
                print("residues_flipped_time trwał: ", residues_flipped_time)

                # 4) Optymalizacje

                self.optimize_with_restraints(model, code, pdb_name)

                self.optimize_full_structure(model, code)

                optimized_end = time.perf_counter()
                optimized_time = optimized_end - residues_flipped_end

                # 5) Zmiana nazwy plików na poprawne (podmienia pliki)
                self.rename_flipped_files(pdb_name)

                # 6) Czyszczenie plików roboczych
                self.cleanup_working_files(pdb_name)


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