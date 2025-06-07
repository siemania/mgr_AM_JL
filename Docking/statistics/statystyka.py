# -*- coding: utf-8 -*-
import os
import sys
import getopt
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as pl
from AutoDockTools.Docking import Docking

# Ścieżka do folderu z plikami wynikowymi
# RESULTS_DIR = "/home/mateusz/IdeaProjects/mgr_AM_JL/Docking/TEST_SKRYPTU/grid_dock_files"
RESULTS_DIR = "presentation_pdb/out"

def parse_dlg_file(file_path):
    """Parsuje plik .dlg, aby znaleźć wartości RMSD i energii (Delta G)."""
    docking = Docking()
    docking.readDlg(file_path)
    id_pdb = os.path.basename(file_path).split(".")[0]  # ID = nazwa pliku bez rozszerzenia

    energy_values = []
    rmsd_values = []

    for cluster in docking.clusterer_dict.values():
        for docking_instance in cluster.clustering_dict.values():
            for docking_pose in docking_instance.data:
                for conformation in docking_pose.data:
                    rmsd_values.append(conformation.refRMS)

                    energy_values.append(conformation.binding_energy)

    return id_pdb, energy_values, rmsd_values


def extract_best_values():
    """Przetwarza pliki .dlg i zwraca słownik {id_pdb: najlepsza wartość} dla wybranego trybu."""
    energy_results = {}
    rmsd_results = {}

    for file in os.listdir(RESULTS_DIR):
        if file.endswith(".dlg"):
            file_path = os.path.join(RESULTS_DIR, file)
            id_pdb, energy_values, rmsd_values = parse_dlg_file(file_path)

            # Wyznaczanie najniższej energii
            if energy_values:
                energy_results[id_pdb] = min(energy_values)

            # Wyznaczanie najniższej wartości RMSD
            if rmsd_values:
                rmsd_results[id_pdb] = min(rmsd_values)

    return energy_results, rmsd_results


def plot_results(energy_data, rmsd_data, exp_energy):
    """Tworzy i zapisuje wykresy Delta G i histogram RMS."""
    # output_dir = "/home/mateusz/IdeaProjects/mgr_AM_JL/Docking/plots"
    output_dir = "plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Tworzenie wykresu Delta G
    if energy_data:
        pl.figure()
        pl.scatter(energy_data.values(), exp_energy)
        pl.xlabel("Wartosci delta G z dokowania")
        pl.ylabel("Eksperymentalne wartosci delta G")
        pl.title("Porownanie wartosci Delta G wzgledem wartosci ekperymentalnych")
        pl.savefig(os.path.join(output_dir, "delta_g_plot.jpg"))
        pl.close()

    # Tworzenie histogramu wartości RMSD
    if rmsd_data:
        pl.figure()
        pl.hist(list(rmsd_data.values()), bins=10, edgecolor='black')
        pl.xlabel("RMSD")
        pl.ylabel("Liczba wynikow")
        pl.title("Histogram RMSD wzgledem liganda eksperymentalnego")
        pl.savefig(os.path.join(output_dir, "rmsd_histogram.jpg"))
        pl.close()


def main(argv):
    """Obsluga argumentów wywołania skryptu."""

    exp_energy = [-7.5, 5, 6, 67, 345, 211]  # Przykładowa wartość eksperymentalna Delta G


    try:
        opts, _ = getopt.getopt(argv, "e:", ["exp_energy="])
    except getopt.GetoptError:
        print("Użycie: extract_energy.py -e <delta_g_eksperymentalne>")
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-e", "--exp_energy"):
            try:
                exp_energy = float(arg)
            except ValueError:
                print("Niepoprawna wartość eksperymentalnej Delta G. Wprowadź liczbę.")
                sys.exit(2)

    energy_results, rmsd_results = extract_best_values()
    print("Wyniki energii:", energy_results)
    print ("Wyniki RMSD:", rmsd_results)

    plots = plot_results(energy_results, rmsd_results, exp_energy)
    print (plots)

if __name__ == "__main__":
    main(sys.argv[1:])
