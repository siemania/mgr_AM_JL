import subprocess

# Ścieżki do programów i skryptów
# Żeby pakiety się zgadzały trzeba odpalać tego Pythona z MGLTools !!!
python = r"E:\Program Files (x86)\MGLTools-1.5.7\python.exe"
dlg = "1hsg_dock.dlg"
receptor = "1hsg_receptor.pdbqt"
file_name = "1hsg"

# Funkcja pomocnicza do uruchamiania poleceń "one by one"
# Pamiętać o podwójnym cudzysłowiu gdy są spacje!
def run_command(command):
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        raise RuntimeError(f"Błąd podczas wykonywania polecenia: {command}")
        
# Utworzenie konformacji z pliku dlg
# run_command(f'"{python}" write_conformations_from_dlg.py -d {dlg}')
# print("Utworzono wszystkie konformacje liganda")

# Utworzenie wszystkich kompleksów z pliku dlg
run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor}')
print("Utworzono wszystkie konformacje kompleksów")

# Utworzenie najlepszego kompleksu ligand-receptor z pliku dlg
run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor} -b -o {file_name}_bestcomplex')
print("Utworzono najlepszą konformację kompleksu")