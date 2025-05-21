import subprocess
import os
import sys

username = os.environ.get('USERNAME') or os.environ.get('USER')
python = f"C:\\Users\\{username}\\OneDrive\\Dyplom_AM_JL\\Skrypty\\Docking\\MGLTools-1.5.7\\python.exe"

def run_command(command):
    process = subprocess.run(command, shell=True)
    if process.returncode != 0:
        raise RuntimeError(f"Błąd podczas wykonywania polecenia: {command}")
        
print(username)
print(python)
print(os.name)
print(sys.platform)

input("~~~~~~~~Naciśnij dowolny przycisk by zakończyć~~~~~~~~")
