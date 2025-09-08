from pathlib import Path
from sys import exit
import argparse

def load_experiment_data(txt):
    """Wczytaj dane z pliku txt do słownika {ID: DeltaG}"""
    data = {}
    with open(txt, "r", encoding="cp1250") as infile:
        header = infile.readline()  # pomiń nagłówek
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                id_val, delta_g = parts[0], parts[3]
                data[id_val.lower()] = delta_g
    return data


if __name__ == "__main__":
    """
    Na podstawie pliku z wartościami energii po obróbce skryptem extract_energy_values.py na plik BioLip.txt
    wyciąga wartości eksperymentalne kcal/mol na podstawie plików wymienionych w --directory
    """

    parser = argparse.ArgumentParser(description="Wyciągnij 1 i 4 kolumnę dla ID obecnych w folderze z plikami pdb")
    parser.add_argument("-f", "--file", required=True, help="Ścieżka do pliku wejściowego (np. ligands_energy_experiment.txt)")
    parser.add_argument("-d", "--directory", required=True, help="Folder z plikami .pdb (np. pdb_files/)")
    parser.add_argument("-o", "--output", required=True, help="Nazwa pliku wyjściowego (np. results.txt)")
    args = parser.parse_args()

    # Wczytaj dane eksperymentalne
    data = load_experiment_data(args.file)

    # Sprawdź folder z pdb
    folder = Path(args.directory)
    if not folder.exists():
        print(f"Error: directory '{folder}' does not exist.")
        exit(1)

    pdb_files = list(folder.glob("*.pdb"))
	
    if len(pdb_files) == 0:
        pdb_files = list(folder.glob("*.dlg"))

    # Przygotuj ścieżkę wyjściową
    output_path = Path(args.output)
    if not output_path.suffix:
        output_path = output_path.with_suffix(".txt")

    # Zapisz wyniki
    with open(output_path, "w", encoding="utf-8") as outfile:
        outfile.write("ID\tDeltaG (kcal/mol)\n")
        for pdb in pdb_files:
            pdb_code = pdb.stem[:4].lower()  # np. '1hi3' z '1hi3.pdb'
            if pdb_code in data:
                outfile.write(f"{pdb_code}\t{data[pdb_code]}\n")

    print(f"Wyniki zapisano do {output_path}")

    # Usage:
    # python extract_experiment_energy.py -f ligands_energy.txt -d results_pdb -o experiment_energy.txt
