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

        # Ważne parametry
        self.env.edat.dynamic_coulomb = True # Uwzględnia oddziaływanie coulombowskie
        self.env.edat.dynamic_sphere = True # Standard dla Modellera zamiast potencjału Lennarda-Jonesa
        self.env.io.hetatm = True # Będzie mogł wczytać atomy liganda

        # STAŁA DIELEKTRYCZNA - coulomb_switch powoduje smoothing stałej relatywnej:
        # oddaje zmienność środowiska dielektrycznego w białku i wokół niego.
        # 1 - Dla symulacji białka — silne oddziaływania elektrostatyczne maleją wraz z odległością
        # 4 - Pomaga "zablokować" ważne mostki solne i wiązania wodorowe w centrum aktywnym w ich optymalnych pozycjach.
        # 80 - Dla symulacji w wodzie (continuum)

        # Wartość standardowa w modellerze
        # self.env.edat.relative_dielectric, self.env.edat.coulomb_switch = 1.0, (6.5, 7.5)

        # Wartości w trakcie testowania
        # self.env.edat.relative_dielectric = 2.0
        # self.env.edat.relative_dielectric = 4.0
        # self.env.edat.relative_dielectric = 8.0
        # self.env.edat.relative_dielectric = 20.0
        # self.env.edat.relative_dielectric = 80.0

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

    def find_active_site(self, protein_model, ligand_pdb_path, cutoff=5.0):
        """
        Identyfikuje reszty centrum aktywnego na podstawie odległości od liganda.
        Wykorzystuje drzewo KDTree (scipy.spatial) do efektywnego wyszukiwania sąsiadujących reszt.
        """
        # 1. Wczytaj atomy LIGANDA z osobnego pliku
        try:
            # Używamy teraz 'model' (klasa) zamiast 'protein_model' (argument)
            ligand_mdl = model(self.env, file=ligand_pdb_path)
        except Exception as e:
            print(f"BŁĄD: Nie można wczytać pliku liganda: {ligand_pdb_path}")
            print(f"Szczegóły: {e}")
            return Selection() # Zwróć pustą selekcję

        ligand_atoms = list(ligand_mdl.atoms)

        if not ligand_atoms:
            print("Nie znaleziono atomów liganda w pliku:", ligand_pdb_path)
            return Selection()

        # Współrzędne atomów liganda
        ligand_coords = np.array([[a.x, a.y, a.z] for a in ligand_atoms])

        # 2. Wczytaj atomy MODELU (białka) - używamy nowej nazwy 'protein_model'
        model_atoms = [a for a in protein_model.atoms if a.residue.name != "HOH"]
        model_coords = np.array([[a.x, a.y, a.z] for a in model_atoms])

        # 3. Zbuduj KDTree dla modelu i znajdź sąsiadów
        tree = cKDTree(model_coords)

        # query_ball_point zwraca listę list indeksów
        nearby_indices_lists = tree.query_ball_point(ligand_coords, cutoff)

        # 4. Zbierz unikalne RESZTY (nie atomy)
        nearby_residues = set()
        for idx_list in nearby_indices_lists:
            for i in idx_list:
                nearby_residues.add(model_atoms[i].residue)

        # 5. Stwórz selekcję ATOMÓW należących do tych reszt
        active_site = Selection(*nearby_residues)

        print(f"Znaleziono {len(nearby_residues)} reszt w centrum aktywnym (w promieniu {cutoff}Å od liganda)")
        return active_site

    def optimize_active_site_sidechains(self, model, pdb_code):
        """
        Optymalizuje tylko atomy łańcuchów bocznych w centrum aktywnym.
        Używa stabilnej metody (CG -> MD -> CG) aby uniknąć "eksplozji".
        Wszystkie pozostałe atomy (szkielet + reszty poza centrum aktywnym) są MROŻONE.
        Wszystko co w centrum aktywnym, pozostaje ruchome.
        """
        print("Rozpoczynanie stabilnej optymalizacji łańcuchów bocznych (CG-MD-CG)...")

        # 1. Identyfikacja centrum aktywnego
        ligand_path = os.path.join(self.project_root, "ligands", f"{pdb_code}_ligand.pdb")
        if not os.path.exists(ligand_path):
            print(f"Ostrzeżenie: Nie znaleziono pliku liganda: {ligand_path}. Pomijanie optymalizacji centrum aktywnego.")
            return

        active_site_residues = self.find_active_site(model, ligand_path, cutoff=7.0)
        if not active_site_residues:
            print("Ostrzeżenie: Nie znaleziono centrum aktywnego. Pomijanie optymalizacji.")
            return

        # 2. Tworzenie selekcji
        backbone_atoms_names = ("N", "CA", "C", "O")

        # active_site_residues to już jest selekcja atomów.
        sidechain_atoms_to_optimize = Selection(*[
            atom for atom in active_site_residues
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

        # 4. Stabilna Optymalizacja (CG -> Annealing MD -> CG)

        # Definicja optymalizatorów
        cg = conjugate_gradients(output='REPORT')
        md = molecular_dynamics(output='REPORT')
        trace_file = open(pdb_code + ".active_site_sidechains_opt.log", "w")

        model.restraints.make(sidechain_atoms_to_optimize, restraint_type='stereo', spline_on_site=False)

        # Krok A: Wstępne CG (usuwa najgorsze zderzenia)
        print("   Krok A: Wstępna minimalizacja CG (30 iteracji)...")
        cg.optimize(sidechain_atoms_to_optimize, max_iterations=30, actions=[actions.trace(10, trace_file)])

        # Krok B: Symulowane Wyżarzanie (Simulated Annealing)
        print("   Krok B: Wyżarzanie MD (600K -> 300K -> 150K -> 50K)...")
        for temp, times in zip([600.0, 300.0, 150.0, 50.0], [150, 75, 75, 75]):
            md.optimize(sidechain_atoms_to_optimize,
                        temperature=temp,
                        max_iterations=times, # Więcej iteracji
                        actions=[actions.trace(10, trace_file)])

        # Krok C: Finałowe CG (dopracowanie geometrii)
        print("   Krok C: Finałowa minimalizacja CG (100 iteracji)...")
        cg.optimize(sidechain_atoms_to_optimize, max_iterations=100, actions=[actions.trace(10, trace_file)])

        trace_file.close()

        # 5. Jawne ODMROŻENIE
        print("Odmrażanie wszystkich zamrożonych atomów.")
        for atom in atoms_to_freeze:
            atom.fix = False

        # # DEBUG
        # energy_result = Selection(model).energy()
        # print(energy_result[0], "Całkowita")
        # print(energy_result[1], "Terms")
        # print(energy_result[1][physical.coulomb], "Coulomb")

        final_path = os.path.join(self.output_path, pdb_code + ".pdb")
        model.write(file=final_path)
        print(f">>> Zoptymalizowano łańcuchy boczne w centrum aktywnym i zapisano: {final_path}")

    def optimize_hydrogens(self, model, pdb_code):
        """
        Przeprowadza stabilną optymalizację (CG -> Annealing MD -> CG) geometrii tylko atomów wodoru.
        Wykorzystuje symulowane wyżarzanie (wysoka temperatura), aby pozwolić terminalnym grupom ("luźnym głowom")
        na swobodną rotację i znalezienie optymalnych pozycji.
        Wszystkie atomy ciężkie są mrożone na czas optymalizacji.
        """
        hydrogen_atoms = Selection(*[atom for atom in model.atoms if atom.name.strip().startswith('H')])
        heavy_atoms = Selection(*[atom for atom in model.atoms if not atom.name.strip().startswith('H')])

        if not hydrogen_atoms:
            print("Ostrzeżenie: Nie znaleziono atomów wodoru do optymalizacji.")
            return

        # 1. ZAMROŻENIE atomów ciężkich (blokuje translacje, ale nie rotację)
        print(f"Mrożenie {len(heavy_atoms)} atomów ciężkich.")
        for atom in heavy_atoms:
            atom.fix = True

        # 2. Stabilna Optymalizacja (CG -> Annealing MD -> CG)
        cg = conjugate_gradients(output='REPORT')
        trace_file = open(pdb_code + ".hydrogens_opt.log", "w")

        model.restraints.make(hydrogen_atoms, restraint_type='stereo', spline_on_site=False)

        print("   Minimalizacja CG wodorów (100 iteracji)...")
        cg.optimize(hydrogen_atoms, max_iterations=100, actions=[actions.trace(10, trace_file)])

        trace_file.close()

        # 3. ODMROŻENIE atomów ciężkich
        print("Odmrażanie atomów ciężkich.")
        for atom in heavy_atoms:
            atom.fix = False

        # # DEBUG
        # energy_result = Selection(model).energy()
        # print(energy_result[0], "Całkowita")
        # print(energy_result[1], "Terms")
        # print(energy_result[1][physical.coulomb], "Coulomb")

        # Nazwa pliku wyjściowego
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

            # Optymalizacja łańcuchów bocznych w centrum aktywnym
            # (na modelu, który ma już zoptymalizowane wodory)
            self.optimize_active_site_sidechains(model, pdb_code)

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