#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, argparse, time, fnmatch
from shutil import copyfile
from tqdm import tqdm
import numpy as np
import flip_res

from modeller import *
from modeller import environ
from modeller.optimizers import molecular_dynamics, actions, conjugate_gradients
from modeller.scripts import complete_pdb
from scipy.spatial import cKDTree


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
        self.input_path = os.path.join(project_root, 'pdb_files') # Folder z wejściowymi plikami PDB
        self.output_path = os.path.join(project_root, 'fixed_pdb') # Folder na naprawione pliki PDB
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
        # Work path delete files
        for file in os.listdir(self.work_path):
            if file.startswith(pdb_code):
                file_path = os.path.join(self.work_path, file)
                try:
                    os.remove(file_path)
                    print(f"Usunięto pliki robocze: {file}")
                except Exception as e:
                    print(f"Error przy usuwaniu {file}: {e}")

        # Fixed_pdb path delete files
        for file in os.listdir(self.output_path): # Usuwa wszystkie pliki z przedrostkami w *_.pdb
            if fnmatch.fnmatch(file, "????_*.pdb"):
                file_path = os.path.join(self.output_path, file)
                try:
                    os.remove(file_path)
                    print(f"Usunięto niepotrzebny plik: {file}")
                except Exception as e:
                    print(f"Error przy usuwaniu {file}: {e}")


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

        print("alignment.ali został zapisany")
        self.aln = aln


    def find_active_site(self, model, ligand_pdb_path, cutoff=5.0):
        """
        Identyfikuje reszty centrum aktywnego na podstawie odległości od liganda.
        Wykorzystuje drzewo KDTree (scipy.spatial) do efektywnego wyszukiwania sąsiadujących reszt.

        Parametry
        ----------
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


    def optimize_heavy_atoms(self, model, pdb_code):
        """
        Przeprowadza lokalną optymalizację geometrii atomów ciężkich (bez H)
        przy użyciu metody gradientów sprzężonych (Conjugate Gradients).
        """
        heavy_atoms = selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])
        model.restraints.make(heavy_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        cg.optimize(heavy_atoms, max_iterations=30)
        final_path = os.path.join(self.output_path, pdb_code + ".pdb")
        model.write(file=final_path)
        print(f">>> Zoptymalizowano atomy ciężkie i zapisano: {final_path}")


    def optimize_with_restraints(self, model, pdb_code):
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
        ligand_path = os.path.join(self.project_root, "ligands", f"{pdb_code}_ligand.pdb")
        active_site = self.find_active_site(model, ligand_path)

        # Więzy: backbone + stereo restraints
        model.restraints.make(backbone, restraint_type='stereo', spline_on_site=False)

        # Zamrożenie centrum aktywnego
        for atom in active_site:
            atom.fix = True

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(pdb_code + ".restrained_opt.log", "w")

        # Gradienty + MD + Gradienty
        cg.optimize(all_atoms, max_iterations=100, actions=[actions.trace(10, trace_file)])
        md.optimize(all_atoms, temperature=300, max_iterations=300,
                    actions=[actions.write_structure(50, pdb_code + ".dyn%04d.pdb"),
                             actions.trace(20, trace_file)])
        cg.optimize(all_atoms, max_iterations=200, actions=[actions.trace(10, trace_file)])

        final_path = os.path.join(self.output_path, pdb_code + "_restrained.pdb")
        model.write(file=final_path)
        print(f">>> Zapisano zoptymalizowaną strukturę (z więzami): {final_path}")


    def optimize_full_structure(self, model, pdb_code):
        """
        Wykonuje pełną optymalizację modelu:
            - relaksacja struktury za pomocą gradientów sprzężonych,
            - symulacja dynamiki molekularnej (MD) z rampowaniem temperatury,
            - ponowna optymalizacja po MD.
        Proces minimalizuje energię układu, poprawiając geometrię lokalną i globalną.
        """
        model.write(file=pdb_code + ".ini")

        all_atoms = Selection(model)
        model.restraints.make(all_atoms, restraint_type='stereo', spline_on_site=False)

        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(pdb_code + ".full_opt.log", "w")

        print("Pełna optymalizacja struktury (CG-MD-CG) ...")

        # Gradienty sprzężone
        cg.optimize(all_atoms, max_iterations=50, actions=[actions.trace(10, trace_file)])

        # Dynamika molekularna z rampowaniem temperatury
        for temp in [150, 300]:
            md.optimize(all_atoms, temperature=temp, max_iterations=100,
                        actions=[actions.write_structure(20, pdb_code + ".dyn%04d.pdb"),
                                 actions.trace(20, trace_file)])

        # Ponowna optymalizacja
        cg.optimize(all_atoms, max_iterations=50, actions=[actions.trace(10, trace_file)])


        final_path = os.path.join(self.output_path, pdb_code + ".pdb")
        model.write(file=final_path)
        print(f">>> Zoptymalizowano pełną strukturę i zapisano: {final_path}")


    def optimize_active_site_sidechains(self, model, pdb_code):
        """
        Optymalizuje tylko atomy łańcuchów bocznych w centrum aktywnym.
        Używa stabilnej metody (CG -> MD -> CG) aby uniknąć "eksplozji".
        Wszystkie pozostałe atomy (szkielet + reszty poza centrum) są MROŻONE.
        """
        print("Rozpoczynanie stabilnej optymalizacji łańcuchów bocznych (CG-MD-CG)...")

        # 1. Identyfikacja centrum aktywnego
        ligand_path = os.path.join(self.project_root, "ligands", f"{pdb_code}_ligand.pdb")
        if not os.path.exists(ligand_path):
            print(f"Ostrzeżenie: Nie znaleziono pliku liganda: {ligand_path}. Pomijanie optymalizacji centrum aktywnego.")
            return

        active_site_residues = self.find_active_site(model, ligand_path)
        if not active_site_residues:
            print("Ostrzeżenie: Nie znaleziono centrum aktywnego. Pomijanie optymalizacji.")
            return

        # 2. Tworzenie selekcji
        backbone_atoms_names = ("N", "CA", "C", "O")
        sidechain_atoms_to_optimize = Selection(*[
            atom for res in active_site_residues.residues
            for atom in res.atoms
            if atom.name.strip() not in backbone_atoms_names
        ])

        if not sidechain_atoms_to_optimize:
            print("Ostrzeżenie: Nie znaleziono atomów łańcuchów bocznych w centrum aktywnym.")
            return

        atoms_to_optimize_set = set(sidechain_atoms_to_optimize)
        atoms_to_freeze = Selection(*[
            atom for atom in model.atoms
            if atom not in atoms_to_optimize_set
        ])

        print(f"Wybrano {len(sidechain_atoms_to_optimize)} atomów łańcuchów bocznych do optymalizacji.")
        print(f"Mrożenie {len(atoms_to_freeze)} pozostałych atomów.")

        # 3. Jawne ZAMROŻENIE
        for atom in atoms_to_freeze:
            atom.fix = True

        # 4. Stabilna Optymalizacja (CG -> MD -> CG)

        # Definicja optymalizatorów
        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(pdb_code + ".active_site_sidechains_opt.log", "w")

        model.restraints.make(sidechain_atoms_to_optimize, restraint_type='stereo', spline_on_site=False)

        # Krok A: Wstępne CG (usuwa najgorsze zderzenia)
        print("   Krok A: Wstępna minimalizacja CG (50 iteracji)...")
        cg.optimize(sidechain_atoms_to_optimize, max_iterations=50, actions=[actions.trace(10, trace_file)])

        # Krok B: Relaksacja MD (pozwala atomom "uciec" ze zderzeń)
        print("   Krok B: Relaksacja MD (100 iteracji @ 300K)...")
        md.optimize(sidechain_atoms_to_optimize,
                    temperature=300.0,
                    max_iterations=100,
                    actions=[actions.trace(20, trace_file)])

        # Krok C: Finałowe CG (dopracowanie geometrii)
        print("   Krok C: Finałowa minimalizacja CG (150 iteracji)...")
        cg.optimize(sidechain_atoms_to_optimize, max_iterations=150, actions=[actions.trace(10, trace_file)])

        trace_file.close()

        # 5. Jawne ODMROŻENIE
        print("Odmrażanie wszystkich zamrożonych atomów.")
        for atom in atoms_to_freeze:
            atom.fix = False

        final_path = os.path.join(self.output_path, pdb_code + ".pdb")
        model.write(file=final_path)
        print(f">>> Zoptymalizowano łańcuchy boczne w centrum aktywnym i zapisano: {final_path}")


    def optimize_hydrogens(self, model, pdb_code):
        """
        Przeprowadza stabilną optymalizację (CG-MD-CG) geometrii tylko atomów wodoru.
        Wszystkie atomy ciężkie są jawnie MROŻONE na czas optymalizacji.
        """
        # Selekcje
        hydrogen_atoms = Selection(*[atom for atom in model.atoms if atom.name.strip().startswith('H')])
        heavy_atoms = Selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])

        if not hydrogen_atoms:
            print("Ostrzeżenie: Nie znaleziono atomów wodoru do optymalizacji.")
            return

        # 1. ZAMROŻENIE atomów ciężkich
        print(f"Mrożenie {len(heavy_atoms)} atomów ciężkich.")
        for atom in heavy_atoms:
            atom.fix = True

        print("Rozpoczynanie optymalizacji atomów wodoru (CG-MD-CG)...")

        # 2. Stabilna Optymalizacja (CG -> MD -> CG)
        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(pdb_code + ".hydrogens_opt.log", "w")

        model.restraints.make(hydrogen_atoms, restraint_type='stereo', spline_on_site=False)

        # Krok A: Wstępne CG
        print("   Krok A: Wstępne CG (30 iteracji)...")
        cg.optimize(hydrogen_atoms, max_iterations=30, actions=[actions.trace(10, trace_file)])

        # Krok B: Relaksacja MD
        print("   Krok B: Relaksacja MD (50 iteracji 300K)...")
        md.optimize(hydrogen_atoms,
                    temperature=300.0,
                    max_iterations=50,
                    actions=[actions.trace(10, trace_file)])

        # Krok C: Finałowe CG
        print("   Krok C: Finałowe CG (100 iteracji)...")
        cg.optimize(hydrogen_atoms, max_iterations=100, actions=[actions.trace(10, trace_file)])

        trace_file.close()

        # 3. ODMROŻENIE atomów ciężkich
        print("Odmrażanie atomów ciężkich.")
        for atom in heavy_atoms:
            atom.fix = False

        final_path = os.path.join(self.output_path, pdb_code + ".pdb")
        model.write(file=final_path)
        print(f">>> Zoptymalizowano atomy wodoru i zapisano: {final_path}")


    def fill_missing_residues_and_atoms(self, files_list=None):
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
        # Get files from the input folder, if provided.
        files = files_list if files_list else next(os.walk(self.input_path))[2]

        for filename in tqdm(files):
            print(f"Przetwarzanie pliku: {filename}")
            start = time.perf_counter()

            pdb_path = os.path.join(self.input_path, filename)
            pdb_code = os.path.splitext(filename)[0]
            temple_name = f"{pdb_code}_fill"

            os.chdir(self.work_path)

            if not os.path.exists(filename):
                print(f"Nie znaleziono {filename}, szukanie w {pdb_path}")
                copyfile(pdb_path, filename)

            self.prepare_alignment(self.env, pdb_code, temple_name)

            prepare_alignment_end = time.perf_counter()
            prepare_alignment_time = prepare_alignment_end - start

            hydrogens_added_end = time.perf_counter()
            hydrogens_added_time = hydrogens_added_end - prepare_alignment_end

            # 1) Wczytanie struktury PDB do uzupełnienia
            # (dodawanie wodorów i ?brakujących atomów?)
            model = complete_pdb(self.env, pdb_path)
            hydrogenated_file = os.path.join(self.work_path, f"{pdb_code}_withH.pdb")
            model.write(file=hydrogenated_file)
            print("Modeller automatycznie dodał wodory (top_allh.lib).")

            model_loaded_end = time.perf_counter()
            model_loaded_time = model_loaded_end - hydrogens_added_end

            # 2) Flipy odpowiednich reszt aminokwasowych
            flipper = flip_res.ResidueFlipper(hydrogenated_file, self.aln)
            flipped_file = flipper.run() # Returns the path to a flipped file.

            # Wczytanie ponownie modelu po flipowaniu
            if os.path.exists(flipped_file):
                model = complete_pdb(self.env, flipped_file)
            else:
                model = complete_pdb(self.env, hydrogenated_file)

            residues_flipped_end = time.perf_counter()
            residues_flipped_time = residues_flipped_end - model_loaded_end

            # 3) Optymalizacje

            # Optymalizacja tylko wodorów
            self.optimize_hydrogens(model, pdb_code)

            # Następnie optymalizujemy łańcuchy boczne w centrum aktywnym
            # (na modelu, który ma już zoptymalizowane wodory)
            # self.optimize_active_site_sidechains(model, pdb_code)

            # self.optimize_with_restraints(model, pdb_code)
            self.optimize_full_structure(model, pdb_code)

            optimized_end = time.perf_counter()
            optimized_time = optimized_end - residues_flipped_end

            # 6) Czyszczenie plików roboczych
            self.cleanup_working_files(pdb_code)

            # 7) Licznik
            files_cleaned_end = time.perf_counter()
            files_cleaned_time = files_cleaned_end - optimized_end

            file_stats = {
                "prepare_alignment_time": prepare_alignment_time,
                "hydrogens_added_time": hydrogens_added_time,
                "model_loaded_time": model_loaded_time,
                "residues_flipped_time": residues_flipped_time,
                "optimized_time": optimized_time,
                "files_cleaned_time": files_cleaned_time
            }

            for key, value in file_stats.items():
                print(f"{key}: {round(value, 2)} s")
            end = time.perf_counter()
            print(f"Całość programu: {round(end - start, 2)} s")



if __name__ == '__main__':
    # Parser argumentów: pozwala uruchomić skrypt dla listy plików lub wszystkich.
    parser = argparse.ArgumentParser(description="Naprawianie struktury PDB (reszty, wodory, flipy)")
    parser.add_argument("-f", "--files", type=str, nargs='+', help="Lista plików .pdb do przetworzenia")
    args = parser.parse_args()
    project_root = os.path.dirname(os.path.abspath(__file__)) # File path to the project root.
    processor = PDBModelOptimization(project_root)

    if args.files:
        files_list = [f if f.endswith('.pdb') else f + '.pdb' for f in args.files]
        processor.fill_missing_residues_and_atoms(files_list=files_list)
    else:
        processor.fill_missing_residues_and_atoms()