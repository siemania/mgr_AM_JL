#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from math import atan, degrees

def parse_sections(filepath):
    """
    Wczytuje dwie sekcje danych z pliku txt.
    Zwraca (names1, x1, y1), (names2, x2, y2), gdzie:
      names = kolumna 1 (nazwa liganda)
      x = kolumna 7 (eksperymentalna energia),
      y = kolumna 4 (energia z dokowania).
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    idx = next(i for i, l in enumerate(lines) if l.startswith("#_Better"))

    def parse_block(block):
        names, xs, ys = [], [], []
        for line in block:
            if line.strip() and not line.startswith("#"):
                parts = line.split()
                try:
                    name = parts[0]
                    y = float(parts[3])
                    x = float(parts[6].replace('|', ''))
                    names.append(name)
                    xs.append(x)
                    ys.append(y)
                except (IndexError, ValueError):
                    continue
        return names, np.array(xs), np.array(ys)

    names1, y1, x1 = parse_block(lines[1:idx])
    names2, y2, x2 = parse_block(lines[idx+2:])
    return (names1, x1, y1), (names2, x2, y2)

def main():
    filepath = "energies_after_docking.txt"
    (names1, x1, y1), (names2, x2, y2) = parse_sections(filepath)

    # regresje
    slope1, intercept1, r1, _, _ = linregress(x1, y1)
    slope2, intercept2, r2, _, _ = linregress(x2, y2)

    # obliczenie odchylenia kąta od 45°
    angle1 = degrees(atan(slope1))
    angle2 = degrees(atan(slope2))
    delta1 = abs(angle1 - 45)
    delta2 = abs(angle2 - 45)

    # zakresy osi
    all_x = np.concatenate([x1, x2])
    all_y = np.concatenate([y1, y2])
    min_val = min(all_x.min(), all_y.min()) - 1
    max_val = max(all_x.max(), all_y.max()) + 1

    # rysowanie
    plt.figure(figsize=(10, 6))
    plt.scatter(x1, y1, label='Standard', marker='x')
    plt.plot(x1, slope1 * x1 + intercept1, 'b--',
             label=f'Regresja Standard (R²={r1**2:.2f}, Δθ={delta1:.2f}°)')

    plt.scatter(x2, y2, label='Better', marker='o')
    plt.plot(x2, slope2 * x2 + intercept2, 'g--',
             label=f'Regresja Better (R²={r2**2:.2f}, Δθ={delta2:.2f}°)')

    # podpisy punktów
    for name, x, y in zip(names1, x1, y1):
        plt.text(x + 0.3, y - 0.4, name, fontsize=8, color='blue')
    for name, x, y in zip(names2, x2, y2):
        plt.text(x - 1, y, name, fontsize=8, color='green')

    plt.title('Porównanie energii wiązania: dokowanie vs dane eksperymentalne')
    plt.ylabel('Energia eksperymentalna [kcal/mol]')
    plt.xlabel('Energia z dokowania [kcal/mol]')
    plt.legend(loc='best')
    plt.grid(True, which='both', linestyle=':', linewidth=0.5)

    # ustawienie tych samych zakresów i kroków osi
    plt.xlim(min_val, max_val)
    plt.ylim(min_val, max_val)
    plt.xticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))
    plt.yticks(np.arange(np.floor(min_val), np.ceil(max_val) + 1, 1))

    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()

    out_png = filepath.rsplit('.',1)[0] + '_correlation.png'
    plt.savefig(out_png, dpi=600)
    print(f"Wykres zapisany do: {out_png}")
    plt.show()

if __name__ == "__main__":
    main()
