#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from math import atan, degrees
import argparse


def load_energy_file(filepath):
    """Wczytuje {ID: DeltaG} z pliku txt (2. kolumna)."""
    data = {}
    with open(filepath, "r", encoding="utf-8") as f:
        header = f.readline()  # pomiń nagłówek
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    data[parts[0]] = float(parts[1])
                except ValueError:
                    continue
    return data


def prepare_data(exp_data, dock_data):
    """Dopasuj wartości (X = docking, Y = experiment) tylko dla wspólnych ID."""
    names, xs, ys = [], [], []
    for lig, dock_val in dock_data.items():
        if lig in exp_data:
            names.append(lig)
            xs.append(dock_val)
            ys.append(exp_data[lig])
    return names, np.array(xs), np.array(ys)


def main():
    parser = argparse.ArgumentParser(
        description="Porównanie energii wiązania: dokowanie vs dane eksperymentalne"
    )
    parser.add_argument("-s", "--standard", required=True, help="Plik TXT z wynikami (np. standard)")
    parser.add_argument("-b", "--fixed", required=True, help="Plik TXT z wynikami (np. fixed)")
    parser.add_argument("-f", "--experiment", required=True, help="Plik TXT z wynikami eksperymentalnymi")
    parser.add_argument("-l", "--labels", required=False, help="Nazwy etykiet", nargs=2, default=("Standard", "Fixed"))
    parser.add_argument("-o", "--output", help="Nazwa pliku wyjściowego (PNG)", default="correlation_plot.png")
    args = parser.parse_args()

    if not args.output.endswith(".png"):
        args.output += ".png"

    # Wczytaj dane
    exp_data = load_energy_file(args.experiment)
    std_data = load_energy_file(args.standard)
    fix_data = load_energy_file(args.fixed)

    # Przygotuj wartości (nazwy, punkty itp.)
    names_std, x_std, y_std = prepare_data(exp_data, std_data)
    names_fix, x_fix, y_fix = prepare_data(exp_data, fix_data)

    # Regresje
    slope_std, intercept_std, r_std, _, _ = linregress(x_std, y_std)
    slope_fix, intercept_fix, r_fix, _, _ = linregress(x_fix, y_fix)

    # Obliczenie odchylenia kąta od 45°
    angle_std = degrees(atan(slope_std))
    angle_fix = degrees(atan(slope_fix))
    delta_std = abs(angle_std - 45)
    delta_fix = abs(angle_fix - 45)

    # Zakresy osi
    all_x = np.concatenate([x_std, x_fix])
    all_y = np.concatenate([y_std, y_fix])
    min_val = min(all_x.min(), all_y.min()) - 1
    max_val = max(all_x.max(), all_y.max()) + 1

    # Rysowanie (s - size)
    plt.figure(figsize=(10, 6))

    plt.scatter(x_std, y_std, label=args.labels[0], marker="x", s=4, color="blue")
    plt.plot(x_std,
             slope_std * x_std + intercept_std,
             "b--",
             label=f"Regresja {args.labels[0]} (R²={r_std**2:.2f}, Δθ={delta_std:.2f}°)")

    plt.scatter(x_fix, y_fix, label=args.labels[1], marker="o", s=4, color="green")
    plt.plot(x_fix,
             slope_fix * x_fix + intercept_fix,
             "g--",
             label=f"Regresja {args.labels[1]} (R²={r_fix**2:.2f}, Δθ={delta_fix:.2f}°)")

    # Podpisy punktów
    # for name, x, y in zip(names_std, x_std, y_std):
    #     plt.text(x + 0.3, y - 0.2, name, fontsize=8, color="blue")
    for name, x, y in zip(names_fix, x_fix, y_fix):
        plt.text(x - 0.4, y + 0.2, name, fontsize=4, color="purple")

    # Podpis wykresu
    plt.title("Porównanie energii wiązania: dokowanie vs dane eksperymentalne")
    plt.ylabel("Energia eksperymentalna [kcal/mol]", fontsize=12)
    plt.xlabel("Energia z dokowania [kcal/mol]", fontsize=12)
    plt.legend(loc="best")
    plt.grid(True, linestyle=":", linewidth=0.5)

    # Parametry wykresu
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.xticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))
    plt.yticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))
    plt.tick_params(axis='both', which='major', labelsize=8) # Wielkość wartości osi

    # Tworzy kwadratowy wykres
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()

    # Zapisywanie wykresu
    plt.savefig(args.output, dpi=600)
    print(f"Wykres zapisany do: {args.output}")
    plt.show()


if __name__ == "__main__":
    main()

# Usage:
# python correlation.py -s wyniki_standard.txt -b wyniki_fixed.txt -f wyniki_eksperymentalne_kcalmol.txt -o wynik_korelacji.png
