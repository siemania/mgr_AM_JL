import subprocess
import os

username = os.environ.get('USERNAME') or os.environ.get('USER')
python = f"C:\\Users\\{username}\\OneDrive\\Dyplom_AM_JL\\Skrypty\\Docking\\MGLTools-1.5.7\\python.exe"

def run_command(command):
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        raise RuntimeError(f"Błąd podczas wykonywania polecenia: {command}")
        
os.makedirs("output_test_folder", exist_ok=True)

run_command(f'"{python}" write_conformations_from_dlg.py -d "1hsg_dock.dlg" -o "output_test_folder/complex"')
