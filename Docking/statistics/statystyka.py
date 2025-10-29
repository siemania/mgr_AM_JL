#!/usr/bin/python
import argparse
import subprocess
from pathlib import Path
import os


def main():
    """
    Skrypt łączący wszystkie inne skrypty wykonujące analizę statystyczną:
    1a. extract_energy_values.py -> wyciąga dane z plików .dlg i tworzy pliki .txt
    1b. !!! Wcześniej należy przerobić plik z bazy BioLip.txt przez binding_energy_reader.py
        !!! oraz extract_experiment_energy.py
    2. rmsd_histogram.py -> tworzy histogram RMSD
    3. energy_correlation.py -> tworzy wykres korelacji energii

    Program przyjmuje 2 foldery oraz nazywa je odpowiednio według uznania.

    Usage:
    python statystyka.py -f ligands_exp.txt -d standard/ fixed/ -n "Zbiór_1" "Zbiór_2"
    """

    # Folder statystyka
    statistics_folder = os.path.dirname(os.path.abspath(Path(__file__)))

    # ArgParser
    parser = argparse.ArgumentParser(description="Pipeline do analizy wyników dokowań")
    parser.add_argument("-d", "--directories", nargs="+", required=True, help="Podaj foldery z plikami .dlg (-d standard fixed)")
    parser.add_argument("-n", "--names", nargs="+", required=True, help="Nazwy odpowiadające folderom (-n Standard Fixed)")
    parser.add_argument("-f", "--file", required=True, help="Plik z wartościami eksperymentalnymi (-f energies_exp.txt)")
    parser.add_argument("-i", "--index", required=False, help="Plik z kodami PDB do uwzględnienia w analizie")
    args = parser.parse_args()

    if len(args.directories) != len(args.names):
        raise ValueError("Liczba nazw (-n) musi odpowiadać liczbie folderów (-d)")

    generated_txt_files = []
    group_name = "_".join(args.names)

    # 1. Uruchom extract_energy_values.py dla każdego folderu
    for folder, name in zip(args.directories, args.names):
        if not Path(folder).exists():
            raise NameError("Brak folderu")
        out_file = f"wyniki_{name}.txt"
        print(f"[INFO] Wyciągam wartości energii z folderu {folder} → {out_file}")
        subprocess.run([
            "python", os.path.join(statistics_folder, "extract_energy_values.py"),
            "-d", folder, "-o", out_file
        ], check=True, shell=True)
        generated_txt_files.append(out_file)

    # 2. Uruchom rmsd_histogram.py na wygenerowanych plikach
    if len(generated_txt_files) >= 2:
        print(f"[INFO] Tworzę histogram RMSD dla {generated_txt_files}")
        cmd = [
            "python", os.path.join(statistics_folder, "rmsd_histogram.py"),
            "-s", generated_txt_files[0],
            "-b", generated_txt_files[1],
            "-l", args.names[0], args.names[1],
            "-o", f"histogram_{group_name}.png"
        ]
        if args.index:
            cmd.extend(["-i", args.index])
        subprocess.run(cmd, check=True, shell=True)

    # 3. Uruchom energy_correlation.py na wygenerowanych plikach + eksperymentalny
    if len(generated_txt_files) >= 2:
        print("[INFO] Tworzę wykres korelacji energii")
        cmd = [
            "python", os.path.join(statistics_folder, "energy_correlation.py"),
            "-s", generated_txt_files[0],
            "-b", generated_txt_files[1],
            "-f", args.file,
            "-l", args.names[0], args.names[1],
            "-o", f"correlation_{group_name}.png"
        ]
        if args.index:
            cmd.extend(["-i", args.index])
        subprocess.run(cmd, check=True, shell=True)

    print("[DONE] Analiza zakończona!")

if __name__ == '__main__':
    main()

# Usage
# python statystyka.py -f ligands_exp.txt -d standard/ fixed/ -n "Standard" "Fixed"
