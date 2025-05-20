#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def extract_rmsd(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    idx = next(i for i, l in enumerate(lines) if l.startswith("#_Better"))

    def parse_block(block):
        """
        Pozwala na podstawie # rozdzielić kolumny
        """
        rmsds = []
        for line in block:
            if line.strip() and not line.startswith("#"):
                parts = line.split()
                try:
                    rmsd = float(parts[2])  # kolumna RMSD
                    rmsds.append(rmsd)
                except (IndexError, ValueError):
                    continue
        return rmsds

    rmsd1 = parse_block(lines[1:idx])
    rmsd2 = parse_block(lines[idx+2:])
    return rmsd1, rmsd2

def plot_rmsd_histogram(rmsd1, rmsd2):
    plt.figure(figsize=(6, 6))
    bins = np.linspace(0, max(max(rmsd1), max(rmsd2)) + 1, 8)

    plt.hist(rmsd1, bins=bins, alpha=0.6, label='Standard', color='skyblue', edgecolor='black')
    plt.hist(rmsd2, bins=bins, alpha=0.6, label='Better', color='lightgreen', edgecolor='black')

    plt.xlabel('RMSD GROMACS [Å]')
    plt.ylabel('Liczba struktur')
    plt.title('Histogram wartości RMSD dla receptorów')
    plt.legend()
    plt.grid(True, linestyle=':')
    plt.tight_layout()
    plt.savefig("histogram_rmsd.png", dpi=600)
    plt.show()

def main():
    # Argumenty systemów Unix
    parser = argparse.ArgumentParser(description="Utwórz histogram na podstawie pliku.")
    parser.add_argument("-f", "--file", help="Podaj ścieżki do plików .txt z RMSD", type=str, nargs='*',
                        default=None)
    args = parser.parse_args()
    
    if args.file:
        filepath = sys.argv[1:]
    else:
        filepath = "energies_after_docking.txt"
    rmsd1, rmsd2 = extract_rmsd(filepath)
    plot_rmsd_histogram(rmsd1, rmsd2)


if __name__ == '__main__':
    main()
