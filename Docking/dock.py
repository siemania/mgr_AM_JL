# -*- coding: utf-8 -*-
import os
import argparse
import glob
import sys
import time
from datetime import date

from config import PYTHON_EXEC, PDBQT_FILES, PDB_FILES, GRID_DOCK_FILES, LIGANDS, MAP_GRID_FILES, OUTPUT_FILES #Sciezki dostepu do odpowiednich plikow
from tqdm import tqdm  # Do pobrania pasek postepu
import subprocess
from multiprocessing import Pool
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
        receptor_pdb = str(PDB_FILES / (file_name + ".pdb"))
        ligand_pdb = str(LIGANDS / (file_name + "_ligand.pdb"))
        receptor_pdbqt = str(PDBQT_FILES / (file_name + "_receptor.pdbqt"))
        ligand_pdbqt = str(PDBQT_FILES / (file_name + "_ligand.pdbqt"))
        gpf = str(GRID_DOCK_FILES / (file_name + "_grid.gpf"))
        glg = str(GRID_DOCK_FILES / (file_name + "_grid.glg"))
        dpf = str(GRID_DOCK_FILES / (file_name + "_dock.dpf"))
        dlg = str(GRID_DOCK_FILES / (file_name + "_dock.dlg"))
        fld = str(MAP_GRID_FILES / (file_name + "_receptor.maps.fld"))
        bestcomplex = str(OUTPUT_FILES / (file_name + "_bestcomplex.pdbqt"))

        # =============================================================================
        #                       start - Wlasciwy program
        # =============================================================================
        # Wykonywanie kolejnych krokow dockingu

        python = str(PYTHON_EXEC)

        subprocess.call([
            python, "prepare_receptor4.py",
            "-r", receptor_pdb,
            "-o", receptor_pdbqt,
            "-A", "bonds_hydrogens",
            "-e", "True",
            "-U", "nphs_lps_waters_nonstdres"
        ])
        if not os.path.exists(receptor_pdbqt):
            raise IOError(
                "Nie utworzono pliku receptor_pdbqt (%s). prepare_receptor4.py się nie powiodlo." % receptor_pdbqt)

        #modify_pdbqt_overwrite(receptor_pdbqt, quiet=True)
        print("Receptor %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, "prepare_ligand4.py",
            "-l", ligand_pdb,
            "-o", ligand_pdbqt
        ])
        print("Ligand %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, "prepare_gpf4.py",
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-y",
            "-o", gpf
        ])
        if not os.path.exists(gpf):
            raise IOError("Nie utworzono pliku GPF dla %s, prawdopodobnie blad w prepare_gpf4.py" % file_name)

        modify_gdpf_overwrite(str(GRID_DOCK_FILES / ("%s_grid.gpf" % file_name)), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Grid %s gotowy" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, "autogrid4",
            "-p", gpf,
            "-l", glg
        ])
        print("Mapy energetyczne %s gotowe" % file_name)
        sys.stdout.flush()

        subprocess.call([
            python, "prepare_dpf42.py",
            "-l", ligand_pdbqt,
            "-r", receptor_pdbqt,
            "-o", dpf
        ])
        modify_gdpf_overwrite(str(GRID_DOCK_FILES / ("%s_dock.dpf" % file_name)), quiet=True)  # Poprawki lokalizacyjne w pliku
        print("Parametry do dokowania %s gotowe" % file_name)
        sys.stdout.flush()

        ## print(f"Rozpoczynanie procesu Autodock-GPU dla {file_name} ...", flush=True)
        ## run_command(f'wsl autodock_gpu_128wi --lfile {ligand_pdbqt} --ffile {fld} --import_dpf {dpf} --resnam {file_name}_lig')
        ## print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)

        print("Rozpoczynanie procesu autodock4.exe dla" + file_name)
        subprocess.call([
            python, "autodock",
            "-p", dpf,
            "-l", dlg
        ])
        print("Dokowanie %s zakonczono pomyslnie" %(file_name))
        sys.stdout.flush()

        # run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor_pdbqt} -o output_files\\{file_name}_bestcomplex -b')
        # print(f"Utworzono najlepszy kompleks ligand-receptor {file_name}", flush=True)

        subprocess.call([
            python, "write_all_complexes.py",
            "-r", receptor_pdbqt,
            "-o", bestcomplex,
            "-b"
        ])
        print("Utworzono najlepszy kompleks ligand-receptor" + file_name)
        sys.stdout.flush()

        # =============================================================================
        #                      koniec - Wlasciwy program
        # =============================================================================
        return u" Zakonczono dokowanie: %s" % (file_name)

    except Exception as e:
        return u" Blad w %s: %s" % (file_name, str(e))

if __name__ == '__main__':

    # Pobranie listy plikow PDB za pomocą parsera argparse
    parser = argparse.ArgumentParser(description="Process docking files.")
    parser.add_argument("-f", "--file", help="Podaj sciezki do plikow PDB lub plik .txt z ID", type=str, nargs='*',
                         default=None)
    args = parser.parse_args()

    if args.file:
         if len(args.file) == 1 and args.file[0].endswith('.txt'): # Sprawdza, czy podano jeden plik .txt
             with open(args.file[0], 'r') as f:
                 ids = [line.split()[0].strip() for line in f] # Otwiera plik i pobiera ID z pierwszej kolumny
             pdb_directory = [os.path.join('pdb_files', "%s.pdb" % id) for id in ids] # Utworzy ścieżki do plików PDB na podstawie ID
         else:
             pdb_directory = args.file # Jeśli nie podano pliku .txt, użyje bezpośrednio podanych ścieżek
    else:
         pdb_directory = glob.glob(os.path.join(str(PDB_FILES), "*.pdb"))


    # Ustawienie liczby rownoleglych procesow
    max_workers = 1

    # Tworzenie katalogow jesli nie istnieja
    for path in [PDBQT_FILES, MAP_GRID_FILES, GRID_DOCK_FILES, OUTPUT_FILES]:
        if not os.path.exists(str(path)):
            os.makedirs(str(path))




    # Pomiar czasu
    date_stamp = date.isoformat(date.today())
    start_time = time.time()
    
    # Uruchamianie przetwarzania w puli procesów + pasek postępu
    #pool = Pool(processes=max_workers)

    results = []

    with open("Docking_log_%s.txt" % date_stamp, "w") as log_file:
        for file in tqdm(pdb_directory, desc="Docking Progress"):
            result = process_docking(file)

            print "\n%s" % result
            log_file.write(result + "\n")
            log_file.flush()
        #tqdm_results = tqdm(pool.imap_unordered(process_docking, pdb_directory), total=len(pdb_directory),
                           #desc="Docking Progress")
        #for file in pdb_directory:
            #process_docking(file)

        #for result in tqdm_results:
            #print("\n%s" % result)
            #log_file.write(result + "\n")
            #log_file.flush()

            elapsed = time.time() - start_time
            estimated_time = format_time(elapsed)

            completed = len(results) + 1
            percent = float(completed) / len(pdb_directory) * 100
            log_file.write("Progress: %d/%d (%.2f%%) | Czas: %s\n" % (completed, len(pdb_directory), percent, estimated_time))
            results.append(result)







    print("Czas zakończenia programu: " + format_time(time.time() - start_time))

    # input("~~~~~~~~Naciśnij dowolny przycisk by zakończyć~~~~~~~~")
    