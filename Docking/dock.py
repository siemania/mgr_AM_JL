# Python 3.12.7
import os
import argparse
import glob
import time
from datetime import date
from tqdm import tqdm # Do pobrania pasek postępu
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from modify_parameters import modify_gdpf_overwrite, modify_pdbqt_overwrite

def format_time(seconds):
    """Obliczanie czasu i sformatowanie > hh:mm:ss"""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def run_command(command): # Pamiętać o podwójnym cudzysłowiu gdy są spacje!
    """Funkcja pomocnicza do uruchamiania poleceń
    stdout=subprocess.DEVNULL > wyłącza widoczność strumienia danych z terminala"""
    process = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL)
    if process.returncode != 0:
        raise RuntimeError(f"Błąd podczas wykonywania polecenia: {command}")


# Funkcja przetwarzająca pojedynczy plik PDB 
def process_docking(pdb_file):
    try:
        # Pobranie nazwy pliku bez rozszerzenia
        file_name = os.path.splitext(os.path.basename(pdb_file))[0]
    
        # Ścieżki do plików
        receptor_pdb = os.path.join("pdb_files", f"{file_name}.pdb")
        ligand_pdb = os.path.join("ligands", f"{file_name}_ligand.pdb")
        receptor_pdbqt = os.path.join("pdbqt_files", f"{file_name}_receptor.pdbqt")
        ligand_pdbqt = os.path.join("pdbqt_files", f"{file_name}_ligand.pdbqt")
        gpf = os.path.join("grid_dock_files", f"{file_name}_grid.gpf")
        glg = os.path.join("grid_dock_files", f"{file_name}_grid.glg")
        dpf = os.path.join("grid_dock_files", f"{file_name}_dock.dpf")
        dlg = os.path.join("grid_dock_files", f"{file_name}_dock.dlg")
    
        # Znalezienie nazwy użytkownika komputera
        mgltools_dir = os.path.join(os.getcwd(), "MGLTools")
        python = os.path.abspath(os.path.join(os.getcwd(), "MGLTools", "python.exe"))
        autogrid = "autogrid4.exe"
        autodock = "autodock4.exe"
    # =============================================================================
    #                       start - Właściwy program
    # =============================================================================    
        # Wykonywanie kolejnych kroków dockingu
        # Tutaj -e -U oraz -A mogą generować błędy i nie wiem jak się ich wyzbyć
        run_command(f'"{python}" prepare_receptor4.py -r {receptor_pdb} -o {receptor_pdbqt} -A "checkhydrogens" -e "True" -U "nphs_lps_waters"')
        modify_pdbqt_overwrite(f"{receptor_pdbqt}", quiet=True)
        print(f"Receptor {file_name} gotowy", flush=True)
        
        run_command(f'"{python}" prepare_ligand4.py -l {ligand_pdb} -o {ligand_pdbqt}')
        print(f"Ligand {file_name} gotowy", flush=True)

        run_command(f'"{python}" prepare_gpf4.py -l {ligand_pdbqt} -r {receptor_pdbqt} -y -o {gpf}')
        modify_gdpf_overwrite(f"grid_dock_files/{file_name}_grid.gpf", quiet=True) # Poprawki lokalizacyjne w pliku
        print(f"Grid {file_name} gotowy", flush=True)

        run_command(f'{autogrid} -p {gpf} -l {glg}')
        print(f"Mapy energetyczne {file_name} gotowe", flush=True)

        run_command(f'"{python}" prepare_dpf42.py -l {ligand_pdbqt} -r {receptor_pdbqt} -o {dpf}')
        modify_gdpf_overwrite(f"grid_dock_files/{file_name}_dock.dpf", quiet=True) # Poprawki lokalizacyjne w pliku
        print(f"Parametry do dokowania {file_name} gotowe", flush=True)

        # print(f"Rozpoczynanie procesu autodock4.exe dla {file_name} ...", flush=True)
        # run_command(f'{autodock} -p {dpf} -l {dlg}')
        # print(f"Dokowanie {file_name} zakonczono pomyslnie", flush=True)
        #
        # run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor_pdbqt} -o output_files\\{file_name}_bestcomplex -b')
        # print(f"Utworzono najlepszy kompleks ligand-receptor {file_name}", flush=True)
    # =============================================================================
    #                      koniec - Właściwy program
    # =============================================================================
        return f"✅ Zakończono: {file_name}"
    
    except Exception as e:
        return f"❌ Błąd w {file_name}: {str(e)}"



if __name__ == '__main__':

    # Pobranie listy plików PDB za pomocą parsera argparse
    parser = argparse.ArgumentParser(description="Process docking files.")
    parser.add_argument("-f", "--file", help="Podaj ścieżki do plików PDB", type=str, nargs='*', default=None)
    args = parser.parse_args()

    if args.file:
        pdb_directory = args.file
    else:
        pdb_directory = glob.glob("pdb_files\\*.pdb")
    
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
    