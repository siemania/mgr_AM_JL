#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg') # Nie tworzy GUI
import matplotlib.pyplot as plt
import argparse


def load_energy_file(filepath, index_pdbs=None):
    """Wczytuje {ID: DeltaG} z pliku txt (2. kolumna)."""
    data = {}
    with open(filepath, "r", encoding="utf-8") as f:
        f.readline()  # pomiń nagłówek
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                try:
                    if index_pdbs is None or parts[0].lower() in index_pdbs:
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
    parser.add_argument("-i", "--index", required=False, help="Plik z kodami PDB do uwzględnienia w analizie")
    args = parser.parse_args()

    if not args.output.endswith(".png"):
        args.output += ".png"

    # Wczytuje kody PDB z pliku (jeśli podano)
    index_pdbs = set()
    if args.index:
        with open(args.index, 'r', encoding="utf-8") as f:
            for line in f:
                index_pdbs.add(line.strip().lower())

    # Wczytaj dane
    exp_data = load_energy_file(args.experiment, index_pdbs if args.index else None)
    std_data = load_energy_file(args.standard, index_pdbs if args.index else None)
    fix_data = load_energy_file(args.fixed, index_pdbs if args.index else None)

    # Przygotuj wartości (nazwy, punkty itp.)
    names_std, x_std, y_std = prepare_data(exp_data, std_data)
    names_fix, x_fix, y_fix = prepare_data(exp_data, fix_data)

    ## Regresje
    slope_std, slope_fix = 1, 1

    # Dla nachylenia = 1 trzeba zmienić b => Wzór: y = 1 * x + b
    intercept_std = np.mean(y_std - x_std)
    intercept_fix = np.mean(y_fix - x_fix)

    # Obliczenie MSE - y: eksperymentalne, x: dokowanie
    mse_std = np.mean((y_std - x_std) ** 2)
    mse_fix = np.mean((y_fix - x_fix) ** 2)
    mse_corr = np.mean((x_std - x_fix) ** 2)

    # Zakresy osi i margines
    all_x = np.concatenate([x_std, x_fix])
    all_y = np.concatenate([y_std, y_fix])
    min_val = min(all_x.min(), all_y.min()) - 1
    max_val = max(all_x.max(), all_y.max()) + 1

    ##############
    ### Wykres ###
    ##############

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10), sharey=True)

    # Linia ukośna
    for ax in (ax1, ax2):
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="grey", linestyle="--", alpha=0.2, linewidth=1)

    # === Panel 1: Standard ===
    ax1.scatter(x_std, y_std, label=args.labels[0], marker="x", s=5, color="blue")
    ax1.plot(x_std,
        slope_std * x_std + intercept_std,
        "b-",
        label=f"Regresja {args.labels[0]}, b={intercept_std:.2f}, MSE={mse_std:.2f})")

    # === Panel 2: Fixed ===
    ax2.scatter(x_fix, y_fix, label=args.labels[1], marker="o", s=5, color="green")
    ax2.plot(x_fix,
        slope_fix * x_fix + intercept_fix,
        "g-",
        label=f"Regresja {args.labels[1]}, b={intercept_fix:.2f}, MSE={mse_fix:.2f}")

    # === Wspólne parametry ===
    for i, ax in enumerate((ax1, ax2)):
        ax.plot([], [], ' ', label=f'MSE między regresjami={mse_corr:.2f}') # Dodaje MSE_corr do legendy
        ax.set_title(f"{args.labels[i]} vs eksperyment")
        ax.set_xlabel("Energia z dokowania [kcal/mol]")
        ax.set_ylabel("Energia eksperymentalna [kcal/mol]")
        ax.grid(True, linestyle=":", linewidth=0.5)
        ax.set_xlim(min_val, max_val)
        ax.set_ylim(min_val, max_val)
        ax.set_xticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))
        ax.set_yticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))
        ax.tick_params(axis='both', which='major', labelsize=8)  # Wielkość wartości osi
        ax.set_aspect("equal", adjustable="box") # Tworzy kwadratowy wykres
        ax.legend(loc='upper left', fontsize=8, framealpha=0.8) # Legenda

    # === Layout ===
    plt.tight_layout(h_pad=1.5)

    # Podpisy punktów (Standardowe, Naprawione)
    for name, x, y in zip(names_std, x_std, y_std):
        if name == '4tw7':
            ax1.text(x + 0.2, y - 0.1, name, fontsize=5, color="blue")
    for name, x, y in zip(names_fix, x_fix, y_fix):
        if name == '4tw7':
            ax2.text(x - 0.2, y + 0.2, name, fontsize=5, color="green")

    # Zapisywanie wykresu
    plt.savefig(args.output, dpi=600)
    print(f"Wykres zapisany do: {args.output}")


if __name__ == "__main__":
    main()

# Usage:
# python energy_correlation.py -s wyniki_standard.txt -b wyniki_fixed.txt -f wyniki_eksperymentalne_kcalmol.txt -o wynik_korelacji.png
