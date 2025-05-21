# Python 3.12.7
import os
import argparse
from glob import glob
import time
from datetime import date
from tqdm import tqdm # Do pobrania pasek postępu
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from modify_parameters import modify_gdpf_overwrite, modify_pdbqt_overwrite, modify_fld_overwrite
from move_files import move_dlg_xml_files

def format_time(seconds):
    """Obliczanie czasu i sformatowanie > hh:mm:ss"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def run_command(command): # Pamiętać o podwójnym cudzysłowiu gdy są spacje!
    """Funkcja pomocnicza do uruchamiania poleceń
    stdout=subprocess.DEVNULL > wyłącza widoczność strumienia danych z terminala"""
    process = subprocess.run(command, stdout=subprocess.DEVNULL, shell=False) # shell=True czasami pomoże/zepsuje coś
    if process.returncode != 0:
        raise RuntimeError(f"Błąd podczas wykonywania polecenia: {command}")


# Funkcja przetwarzająca pojedynczy plik PDB 
def process_docking(pdb_file, commands=None):
    try:
        # Pobranie nazwy pliku bez rozszerzenia
        file_name = os.path.splitext(os.path.basename(pdb_file))[0]
    
        # Ścieżki do plików - absolute path nie zawsze działa, bo tworzy w parametrach pełne ścieżki
        receptor_pdb = os.path.join("pdb_files", f"{file_name}.pdb")
        ligand_pdb = os.path.join("ligands", f"{file_name}_ligand.pdb")
        receptor_pdbqt = os.path.join("pdbqt_files", f"{file_name}_receptor.pdbqt")
        ligand_pdbqt = os.path.join("pdbqt_files", f"{file_name}_ligand.pdbqt")
        gpf = os.path.join("grid_dock_files", f"{file_name}_grid.gpf")
        glg = os.path.join("grid_dock_files", f"{file_name}_grid.glg")
        dpf = os.path.join("grid_dock_files", f"{file_name}_dock.dpf")
        dlg = os.path.join("grid_dock_files", f"{file_name}_ligand.dlg")
        fld = os.path.join("map_grid_files", f"{file_name}_receptor.maps.fld")
    
        # Ścieżki do innych programów (kompatybilność)
        user = os.environ.get("USER") or os.environ.get("USERNAME")
        mgltools_dir = os.path.abspath(os.path.join(os.getcwd(), "MGLTools"))
        python = os.path.abspath(os.path.join(os.getcwd(), "MGLTools", "python.exe"))
        autogrid = os.path.abspath(os.path.join(os.getcwd(), "autogrid4.exe"))
        autodock_gpu_path = glob(f"/home/{user}/**/autodock_*wi", recursive=True)[0]
        autodock = os.path.abspath(os.path.join(os.getcwd(), "autodock4.exe"))
    # =============================================================================
    #                       start - Właściwy program
    # =============================================================================    
        if not commands or 'receptor' in commands:
            # Przygotowanie receptorów
            run_command([
                            python,
                            "prepare_receptor4.py",
                            "-r", receptor_pdb,
                            "-o", receptor_pdbqt,
                            "-A", "checkhydrogens",
                            "-e", "True",
                            "-U", "nphs_lps_waters_nonstdres"
                        ])
            modify_pdbqt_overwrite(receptor_pdbqt, quiet=True)
            print(f"\nReceptor {file_name} gotowy", flush=True)

        if not commands or 'ligand' in commands:
            # Przygotowanie ligandów
            run_command([
                            python,
                            "prepare_ligand4.py",
                            "-l", ligand_pdb,
                            "-o", ligand_pdbqt
                        ])
            print(f"Ligand {file_name} gotowy", flush=True)

        if not commands or 'grid' in commands:
            # Przygotowanie parametrów do grida
            run_command([
                            python,
                            "prepare_gpf4.py",
                            "-l", ligand_pdbqt,
                            "-r", receptor_pdbqt,
                            "-y",
                            "-o", gpf
                        ])
            modify_gdpf_overwrite(f"grid_dock_files/{file_name}_grid.gpf", quiet=True)  # Poprawki lokalizacyjne w pliku
            print(f"Parametry do grida {file_name} gotowe", flush=True)

        if not commands or 'autogrid' in commands:
            # Przygotowanie map energetycznych
            run_command([
                            autogrid,
                            "-p", gpf,
                            "-l", glg
                        ])
            modify_fld_overwrite(fld, quiet=True)
            print(f"Mapy energetyczne {file_name} gotowe", flush=True)

        if not commands or 'dock' in commands:
            # Przygotowanie parametrów do dokowania
            run_command([
                            python,
                            "prepare_dpf42.py",
                            "-l", ligand_pdbqt,
                            "-r", receptor_pdbqt,
                            "-o", dpf
                        ])
            modify_gdpf_overwrite(f"grid_dock_files/{file_name}_dock.dpf", quiet=True) # Poprawki lokalizacyjne w pliku
            print(f"Parametry do dokowania {file_name} gotowe", flush=True)

        if not commands or 'autodock' in commands:
            print(f"Rozpoczynanie procesu Autodock-GPU dla {file_name} ...", flush=True)
            run_command([
                            # Ścieżka do programu w systemie WSL (np. /home/user/AutoDock-GPU/...)
                            autodock_gpu_path,
                            "--lfile", ligand_pdbqt,
                            "--ffile", fld,
                            "--nrun", "10",
                        ])
            move_dlg_xml_files("pdbqt_files", "grid_dock_files")
            print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)

        # DLA KOMPATYBILNOŚCI
        # if not commands or 'autodock' in commands:
        #     print(f"Rozpoczynanie procesu autodock4.exe dla {file_name} ...", flush=True)
        #     run_command([
        #                     autodock,
        #                     "-p", dpf,
        #                     "-l", dlg
        #                 ])
        #
        #     print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)

        if not commands or 'complex' in commands:
            run_command([
                            python,
                            "write_all_complexes.py",
                            "-d", dlg,
                            "-r", receptor_pdbqt,
                            "-o", f"output_files/{file_name}_bestcomplex",
                            "-b"
                        ])
            print(f"Utworzono najlepszy kompleks ligand-receptor {file_name}", flush=True)
    # =============================================================================
    #                      koniec - Właściwy program
    # =============================================================================
        return f"✅ Zakończono: {file_name}"
    
    except Exception as e:
        return f"❌ Błąd w {file_name}: {str(e)}"


if __name__ == '__main__':
    # Pobranie listy plików PDB za pomocą parsera argparse
    parser = argparse.ArgumentParser(description="Process docking files.")
    parser.add_argument("-f", "--file", help="Podaj ścieżki do plików PDB lub plik .txt z ID", type=str, nargs='*',
                        default=None)
    parser.add_argument("-s", "--select_command",
                        help="Wybierz komendy do wykonania: receptor, ligand, grid, autogrid, dock, autodock, complex.",
                        type=str, nargs="+", choices=['receptor', 'ligand', 'grid', 'autogrid', 'dock', 'autodock', 'complex'])

    args = parser.parse_args()

    # Pozwala wykorzystać 1 kolumnę z pliku .txt lub wypisać własne ścieżki do plików lub wykorzysta folder pdb_files
    if args.file:
        if len(args.file) == 1 and args.file[0].endswith('.txt'): # Sprawdza, czy podano jeden plik .txt
            with open(args.file[0], 'r') as f:
                ids = [line.split()[0].strip() for line in f] # Otwiera plik i pobiera ID z pierwszej kolumny
            pdb_directory = [os.path.join('pdb_files', f'{id}.pdb') for id in ids] # Utworzy ścieżki do plików PDB na podstawie ID
        else:
            pdb_directory = args.file # Jeśli nie podano pliku .txt, użyje bezpośrednio podanych ścieżek
    else:
        pdb_directory = glob("pdb_files\\*.pdb")
    
    # Ustawienie liczby równoległych procesów
    max_workers = 4
    
    # Tworzenie katalogów jeśli nie istnieją
    os.makedirs("pdbqt_files", exist_ok=True)
    os.makedirs("map_grid_files", exist_ok=True)
    os.makedirs("grid_dock_files", exist_ok=True)
    os.makedirs("output_files", exist_ok=True)
    
    # Pomiar czasu
    date_stamp = date.isoformat(date.today())
    start_time = time.time()
    
    # Uruchamianie przetwarzania w puli procesów + pasek postępu
    with ProcessPoolExecutor(max_workers=max_workers) as executor, tqdm(total=len(pdb_directory), desc="Docking Progress") as progress:
        future_to_file = {executor.submit(process_docking, pdb_file, args.select_command): pdb_file for pdb_file in
                          pdb_directory}

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
    