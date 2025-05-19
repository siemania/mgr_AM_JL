# -*- coding: utf-8 -*-
import os
import glob
import sys
import time
from __builtin__ import raw_input
from datetime import date

from config import PYTHON_EXEC, PDBQT_FILES, PDB_FILES, GRID_DOCK_FILES, LIGANDS, MAP_GRID_FILES, OUTPUT_FILES
from tqdm import tqdm  # Do pobrania pasek postepu
import subprocess


# from modify_parameters import modify_gdpf_overwrite


def format_time(seconds):
    """Obliczanie czasu i sformatowanie > hh:mm:ss"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "%02d:%02d:%02d" % (int(hours), int(minutes), int(seconds))


# Funkcja przetwarzajaca pojedynczy plik PDB
def process_docking(pdb_file):
    try:
        # Pobranie nazwy pliku bez rozszerzenia
        file_name = os.path.splitext(os.path.basename(pdb_file))[0]

        # Sciezki do plikow
        receptor_pdb = str(PDB_FILES / (file_name + ".pdb"))
        ligand_pdb = str(LIGANDS / (file_name + "_ligand.pdb"))
        receptor_pdbqt = str(PDBQT_FILES / (file_name + "_receptor.pdbqt"))
        ligand_pdbqt = str(PDBQT_FILES / (file_name + "_ligand.pdbqt"))
        gpf = str(GRID_DOCK_FILES / (file_name + "_grid.gpf"))
        glg = str(GRID_DOCK_FILES / (file_name + "_grid.glg"))
        dpf = str(GRID_DOCK_FILES / (file_name + "_dock.dpf"))
        dlg = str(GRID_DOCK_FILES / (file_name + "_dock.dlg"))

        # =============================================================================
        #                       start - Wlasciwy program
        # =============================================================================
        # Wykonywanie kolejnych krokow dockingu

        python = str(PYTHON_EXEC)

        subprocess.call([
            python, str(TEST_SKRYPTU_DIR / "prepare_receptor4.py"),
            "-r", receptor_pdb,
            "-o", receptor_pdbqt,
            "-A", "bonds_hydrogens",
            "-e", "True",
            "-U", "nphs","waters","lps"
        ])

        print("Receptor %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, str(TEST_SKRYPTU_DIR / "prepare_ligand4.py"),
            "-l", ligand_pdb,
            "-o", ligand_pdbqt
        ])
        print("Ligand %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, str(TEST_SKRYPTU_DIR / "prepare_gpf4.py"),
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-o", gpf
        ])
        # modify_gdpf_overwrite(os.path.join(GRID_DOCK_FILES, "%s_grid.gpf" % file_name), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Grid %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            "autogrid4",
            "-p", gpf,
            "-l", glg
        ])
        print("Mapy energetyczne %s gotowe" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, str(TEST_SKRYPTU_DIR / "prepare_dpf42.py"),
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-o", dpf
        ])
        # modify_gdpf_overwrite(os.path.join(GRID_DOCK_FILES, "%s_dock.dpf" % file_name), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Parametry do dokowania %s gotowe" % file_name)
        sys.stdout.flush()

        # print("Rozpoczynanie procesu autodock4.exe dla %s ..." % file_name)
        # sys.stdout.flush()
        # run_command('%s -p %s -l %s' % (autodock, dpf, dlg))
        # print("Dokowanie %s zakonczono pomyslnie" % file_name)
        # sys.stdout.flush()

        # run_command('"%s" write_all_complexes.py -d %s -r %s -o output_files\\%s_bestcomplex -b' % (python, dlg, receptor_pdbqt, file_name))
        # print("Utworzono najlepszy kompleks ligand-receptor %s" % file_name)
        # sys.stdout.flush()
        # =============================================================================
        #                      koniec - Wlasciwy program
        # =============================================================================
        return u" Zakonczono: %s" % (file_name)

    except Exception as e:
        return u" Blad w %s: %s" % (file_name, str(e))


if __name__ == '__main__':

    # Pobranie listy plikow PDB
    pdb_directory = glob.glob(str(PDB_FILES / "*.pdb"))

    # Ustawienie liczby rownoleglych procesow
    max_workers = 4

    # Tworzenie katalogow jesli nie istnieja
    for path in [PDBQT_FILES, MAP_GRID_FILES, GRID_DOCK_FILES, OUTPUT_FILES]:
        if not os.path.exists(str(path)):
            os.makedirs(str(path))

    # Pomiar czasu
    date_stamp = date.isoformat(date.today())
    start_time = time.time()

    # Uruchamianie przetwarzania w puli procesow + pasek postepu
    #  pool = Pool(processes=max_workers)
    results = []
    log_path = "Docking_log_%s.txt" % date_stamp
    log_file = open(log_path, "w")

    progress = tqdm(total=len(pdb_directory), desc="Docking Progress")

    for file in  pdb_directory:
        # print("\n" + result)
        # log_file.write(result + "\n")
        # log_file.flush()
        # progress.update(1)

        process_docking(file)

        elapsed = progress.format_dict.get("elapsed", 0)
        estimated_time = format_time(elapsed)
        log_file.write("Progress: %d/%d (%.2f%%) | Czas: %s\n" %
                       (progress.n, progress.total, (float(progress.n) / progress.total) * 100, estimated_time))

    progress.close()
    log_file.close()
    print("Czas zakonczenia programu: %s" % format_time(time.time() - start_time))
    raw_input("~~~~~~~~Nacisnij dowolny przycisk by zakonczyc~~~~~~~~")
