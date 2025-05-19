# -*- coding: utf-8 -*-
import os
import argparse
import glob
import sys
import time
from __builtin__ import raw_input
from datetime import date

from config import PYTHON_EXEC, PDBQT_FILES, PDB_FILES, GRID_DOCK_FILES, LIGANDS, MAP_GRID_FILES, OUTPUT_FILES #Sciezki dostepu do odpowiednich plikow
from tqdm import tqdm  # Do pobrania pasek postepu
import subprocess
from multiprocessing import Pool  # Zamiast `ProcessPoolExecutor`
from modify_parameters import modify_gdpf_overwrite, modify_pdbqt_overwrite

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
        receptor_pdb = os.path.join(PDB_FILES, file_name + ".pdb")
        ligand_pdb = os.path.join(LIGANDS, file_name + "_ligand.pdb")
        receptor_pdbqt = os.path.join(PDBQT_FILES, file_name + "_receptor.pdbqt")
        ligand_pdbqt = os.path.join(PDBQT_FILES, file_name + "_ligand.pdbqt")
        gpf = os.path.join(GRID_DOCK_FILES, file_name + "_grid.gpf")
        glg = os.path.join(GRID_DOCK_FILES, file_name + "_grid.glg")
        dpf = os.path.join(GRID_DOCK_FILES, file_name + "_dock.dpf")
        dlg = os.path.join(GRID_DOCK_FILES, file_name + "_dock.dlg")
        fld = os.path.join(MAP_GRID_FILES, file_name + "_receptor.maps.fld")

        # =============================================================================
        #                       start - Wlasciwy program
        # =============================================================================
        # Wykonywanie kolejnych krokow dockingu
        # Tutaj -e -U oraz -A mogą generowaC bLEdy i nie wiem jak się ich wyzbyC

        subprocess.call([
            PYTHON_EXEC, "prepare_receptor4.py",
            "-r", receptor_pdb,
            "-o", receptor_pdbqt,
            "-A", "checkhydrogens",
            "-e", "True",
            "-U", "nphs_waters_lps_nonstdres"
        ])
        modify_pdbqt_overwrite("%s" % receptor_pdbqt, quiet=True)
        print("Receptor %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            PYTHON_EXEC, "prepare_ligand4.py",
            "-l", ligand_pdb,
            "-o", ligand_pdbqt
        ])
        print("Ligand %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            PYTHON_EXEC, "prepare_gpf4.py",
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-y",
            "-o", gpf
        ])
        modify_gdpf_overwrite(os.path.join(GRID_DOCK_FILES, "%s_grid.gpf" % file_name), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Grid %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            PYTHON_EXEC, "autogrid4",
            "-p", gpf,
            "-l", glg
        ])
        print("Mapy energetyczne %s gotowe" % file_name)
        sys.stdout.flush()

        subprocess.call([
            PYTHON_EXEC, "prepare_dpf42.py",
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-o", dpf
        ])
        modify_gdpf_overwrite(os.path.join(GRID_DOCK_FILES, "%s_dock.dpf" % file_name), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Parametry do dokowania %s gotowe" % file_name)
        sys.stdout.flush()

        ## print(f"Rozpoczynanie procesu Autodock-GPU dla {file_name} ...", flush=True)
        ## run_command(f'wsl autodock_gpu_128wi --lfile {ligand_pdbqt} --ffile {fld} --import_dpf {dpf} --resnam {file_name}_lig')
        ## print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)

        # print(f"Rozpoczynanie procesu autodock4.exe dla {file_name} ...", flush=True)
        # run_command(f'{autodock} -p {dpf} -l {dlg}')
        # print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)
        #
        # run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor_pdbqt} -o output_files\\{file_name}_bestcomplex -b')
        # print(f"Utworzono najlepszy kompleks ligand-receptor {file_name}", flush=True)

        # =============================================================================
        #                      koniec - Wlasciwy program
        # =============================================================================
        return u" Zakonczono: %s" % (file_name)

    except Exception as e:
        return u" Blad w %s: %s" % (file_name, str(e))

if __name__ == '__main__':

    # Pobranie listy plikow PDB za pomocą parsera argparse
    parser = argparse.ArgumentParser(description="Process docking files.")
    parser.add_argument("-f", "--file", help="Podaj ścieżki do plików PDB lub plik .txt z ID", type=str, nargs='*',
                        default=None)
    args = parser.parse_args()

    if args.file:
        if len(args.file) == 1 and args.file[0].endswith('.txt'): # Sprawdza, czy podano jeden plik .txt
            with open(args.file[0], 'r') as f:
                ids = [line.split()[0].strip() for line in f] # Otwiera plik i pobiera ID z pierwszej kolumny
            pdb_directory = [os.path.join('pdb_files', f'{id}.pdb') for id in ids] # Utworzy ścieżki do plików PDB na podstawie ID
        else:
            pdb_directory = args.file # Jeśli nie podano pliku .txt, użyje bezpośrednio podanych ścieżek
    else:
        pdb_directory = glob.glob(os.path.join(PDB_FILES, "*.pdb"))

    # Ustawienie liczby rownoleglych procesow
    max_workers = 4

    # Tworzenie katalogow jesli nie istnieja
    if not os.path.exists(PDBQT_FILES):
        os.makedirs(PDBQT_FILES)
    if not os.path.exists(MAP_GRID_FILES):
        os.makedirs(MAP_GRID_FILES)
    if not os.path.exists(GRID_DOCK_FILES):
        os.makedirs(GRID_DOCK_FILES)
    if not os.path.exists(OUTPUT_FILES):
        os.makedirs(OUTPUT_FILES)

    # Pomiar czasu
    date_stamp = date.isoformat(date.today())
    start_time = time.time()
    
    # Uruchamianie przetwarzania w puli procesów + pasek postępu
    with ProcessPoolExecutor(max_workers=max_workers) as executor, tqdm(total=len(pdb_directory), desc="Docking Progress") as progress:
        future_to_file = {executor.submit(process_docking, pdb_file): pdb_file for pdb_file in pdb_directory}
    
        # Zapis do Docking_log.txt
        with open(f"Docking_log_{date_stamp}.txt", "w", encoding="utf-8") as log_file:
    
            # Monitorowanie zakończonych zadań w czasie rzeczywistym
            for future in as_completed(future_to_file):
                result = future.result()
                print(f"\n{result}") # Wypisuje komunikat o zakończonym pliku
                
                # Zapis do logu
                log_file.write(result + "\n")
                log_file.flush() # Zapis na bieżąco

                # Aktualizacja paska postępu i zapis do pliku
                progress.update(1)
                estimated_time = format_time(progress.format_dict.get("elapsed")) # Pozostały czas bez ms
                log_file.write(
                    f"Progress: {progress.n}/{progress.total} ({(progress.n / progress.total)*100:.2f}%) | Czas: {estimated_time}\n")

    print(f'Czas zakończenia programu: {format_time(time.time() - start_time)}')
    # input("~~~~~~~~Naciśnij dowolny przycisk by zakończyć~~~~~~~~")
    