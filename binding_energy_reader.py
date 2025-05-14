import os
import re
import math
import pandas as pd
import argparse

def compute_delta_g(K_molar, temperature=298.15):
    """
    Oblicza zmianę energii wiązania ΔG w kcal/mol: ΔG = R * T * ln(K)
    R = 0.0019872036 kcal/(mol*K)
    """
    R_kcal = 0.0019872036
    return R_kcal * temperature * math.log(K_molar)


def parse_binding_constant(line):
    """
    Parsuje linię pod kątem Ki, Kd, IC50 lub -log(Ki/Kd).
    Zwraca wartość stałej (M) oraz typ (str).
    """
    match = re.search(r"\b(Ki|Kd|IC50)\s*=\s*([0-9\.Ee+-]+)\s*([munp]?M)\b", line)
    if match:
        typ, val, unit = match.groups()
        K = float(val)
        factors = {'M':1, 'mM':1e-3, 'uM':1e-6, 'nM':1e-9, 'pM':1e-12}
        return K * factors[unit], typ
    match2 = re.search(r"-log\s*\((?:Ki|Kd)\)\s*=\s*([0-9\.Ee+-]+)\b", line)
    if match2:
        pval = float(match2.group(1))
        K = 10 ** (-pval)
        return K, 'pKi'
    return None, None


def process_ligands(csv_path, biolip_txt, output_path, temperature=298.15):
    df = pd.read_csv(csv_path)
    if 'id' not in df.columns:
        raise KeyError("Brak kolumny 'id' w CSV")

    with open(biolip_txt, 'r') as f:
        lines = f.readlines()

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    results = []
    for ligand_id in df['id']:
        pattern = re.compile(rf"^{re.escape(str(ligand_id))}\b")
        line = next((l for l in lines if pattern.search(l)), None)
        if not line:
            msg = "Brak wpisu w BioLiP"
            print(f"{ligand_id}: {msg}", flush=True)
            results.append((ligand_id, None, None, msg))
            continue
        K_molar, typ = parse_binding_constant(line)
        if K_molar is None:
            msg = "Brak stałej wiązania"
            print(f"{ligand_id}: {msg}", flush=True)
            results.append((ligand_id, None, None, msg))
            continue
        delta_g = compute_delta_g(K_molar, temperature)
        results.append((ligand_id, typ, K_molar, f"{delta_g:.3f}"))

    with open(output_path, 'w') as out:
        header = "ID\tTyp\tK (M)\tDeltaG (kcal/mol)\n"
        out.write(header)
        for lig, typ, K, dg in results:
            line = f"{lig}\t{typ or '-'}\t{K if K else '-'}\t{dg if dg else '-'}\n"
            out.write(line)
            # Informacja o zapisie każdego liganda
            print(f"Zapisano dla {lig}: Typ={typ or '-'}, K={K if K else '-'}, DeltaG={dg if dg else '-'}", flush=True)
    print(f"\nGotowe! Wszystkie wyniki zapisane w {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Obliczanie ΔG dla ligandów z BioLiP")
    parser.add_argument('-c', '--csv', default='baza_ids.csv', help='Ścieżka do pliku CSV')
    parser.add_argument('-b', '--biolip', default='BioLiP.txt', help='Ścieżka do pliku BioLiP')
    parser.add_argument('-o', '--output', default=os.path.join('pdb_energy', 'pdb_energy.txt'),
                        help='Plik wynikowy')
    parser.add_argument('-t', '--temp', type=float, default=298.15, help='Temperatura w K')
    args = parser.parse_args()

    process_ligands(args.csv, args.biolip, args.output, args.temp)

if __name__ == '__main__':
    main()

    # Użycie:
    # python binding_energy_reader.py -c baza_ids.csv -b BioLiP.txt -o pdb_energy/ligands_energy.txt