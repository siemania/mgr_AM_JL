#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse

def extract_rmsd(filepath):
    """Wczytuje trzecią kolumnę RMSD z pliku TSV"""
    rmsds = []
    with open(filepath, 'r', encoding="utf-8") as f:
        header = f.readline()  # pomiń nagłówek
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                try:
                    rmsd = float(parts[2])
                    rmsds.append(rmsd)
                except ValueError:
                    continue
    return rmsds


# Standardowe opcje wyboru zapewnia argparse, ale tylko dla tego pliku rozruchowego!
def plot_rmsd_histogram(rmsd1, rmsd2, labels=("Standard", "Fixed"), output="histogram_rmsd.png"):
    plt.figure(figsize=(6, 6))
    bins = np.linspace(0, max(max(rmsd1), max(rmsd2)) + 1, num=12) # num - bins number /4

    plt.hist(rmsd1, bins=bins, alpha=0.6, label=labels[0],
             color='skyblue', edgecolor='black')
    plt.hist(rmsd2, bins=bins, alpha=0.6, label=labels[1],
             color='lightgreen', edgecolor='black')

    plt.xlabel('RMSD GROMACS [Å]')
    plt.ylabel('Liczba struktur')
    plt.title('Histogram wartości RMSD dla receptorów')
    plt.legend()
    plt.grid(True, linestyle=':')
    plt.tight_layout()
    plt.savefig(output, dpi=600)
    print(f"Wykres RMSD zapisany do: {output}")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Utwórz histogram na podstawie dwóch plików TXT z RMSD.")
    parser.add_argument("-s", "--standard", required=True, help="Plik TXT z wynikami standard")
    parser.add_argument("-b", "--fixed", required=True, help="Plik TXT z wynikami fixed")
    parser.add_argument("-l", "--labels", required=False, help="Nazwy etykiet", nargs=2, default=("Standard", "Fixed"))
    parser.add_argument("-o", "--output", required=False, help="Plik wynikowy .png", default="histogram_rmsd.png", type=str)
    args = parser.parse_args()

    # Suffix
    if not args.output.endswith(".png"):
        args.output += ".png"

    rmsd_standard = extract_rmsd(args.standard)
    rmsd_fixed = extract_rmsd(args.fixed)

    plot_rmsd_histogram(rmsd_standard, rmsd_fixed,
                        labels=args.labels, output=args.output)

# Usage:
# python rmsd_hist.py -s wyniki_standard.txt -b wyniki_fixed.txt
