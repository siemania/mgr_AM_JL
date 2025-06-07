from pathlib import Path
import argparse

# Funkcja do wyciągania RMSD i energii z najlepiej zadokowanej konformacji pliku .dlg
def extract_best_conformation_info(dlg_file_path) -> tuple:
    min_energy = float('inf')
    min_rmsd = None

    with open(dlg_file_path, 'r') as file:
        for line in file:
            # Rozpoczęcie bloku konformacji
            if line.endswith("RANKING\n"):
                min_energy = float(line[22:33].strip())
                min_rmsd = float(line[44:61].strip())
                
                return (min_energy, min_rmsd)
            
        raise Exception("Minimum energy and rmsd not found")


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Extract minimum energy and RMSD from .dlg docking files.")
    parser.add_argument("-f", "--file", help="One or more .dlg files to process", type=str, nargs='*')
    parser.add_argument("-d", "--directory", help="Directory containing .dlg files", type=str)
    parser.add_argument("-o", "--output", help="Output file name (with or without .txt extension)", type=str, default="energy_rmsd.txt")
    args = parser.parse_args()

    dlg_files = []

    if args.file:
        dlg_files.extend([Path(f) for f in args.file])
    elif args.directory:
        dlg_folder = Path(args.directory)
        if not dlg_folder.exists():
            print(f"Error: directory '{dlg_folder}' does not exists.")
            exit(1)
        dlg_files.extend(dlg_folder.glob("*.dlg"))
    else:
        # Domyślny folder, jeśli nie podano niczego
        dlg_folder = Path("better_minimize")
        dlg_files.extend(dlg_folder.glob("*.dlg"))

    # Popraw nazwę pliku wyjściowego jeśli trzeba
    output_path = Path(args.output)
    
    if not output_path.suffix:
        output_path = output_path.with_suffix(".txt")

    # Lista wyników
    results = ["PDB\tENERGY\tRMSD"]

    for dlg_file in dlg_files:
        try:
            min_energy, min_rmsd = extract_best_conformation_info(dlg_file)
            pdb_code = dlg_file.stem[:4]
            results.append(f"{pdb_code}\t{min_energy:.2f}\t{min_rmsd:.2f}")
        except Exception as e:
            print(f"Error in file {dlg_file.name}: {e}")

    # Zapisz do pliku
    with open(output_path, "w") as f:
        f.write("\n".join(results))

    print(f"\nData saved in: {output_path.resolve()}")
    
    # Użycie:
    # python extract_docking_info.py -f example.dlg
    # python extract_docking_info.py -f 1abc.dlg 2def.dlg 3ghi.dlg
    # python extract_docking_info.py -d docking_results/
    # python extract_docking_info.py -d docking_results/ -o wyniki_dokowania.txt
