#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Skrypt do minimalizacji receptora przy użyciu OpenBabel.
Wczytuje receptor (PDB, czysty – bez wodorów), dodaje brakujące atomy wodoru,
a następnie przeprowadza dwustopniową minimalizację energii:
  Krok 1: Minimalizacja wyłącznie atomów wodoru (ciężkie atomy są zablokowane).
  Krok 2: Krótka, łagodna minimalizacja całego układu (po usunięciu restrykcji).
Wynik zapisuje do nowego pliku PDB.
"""

import os
from openbabel import openbabel

class ReceptorMinimizerOpenBabel:
    def __init__(self, pdb_file):
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"Plik {pdb_file} nie został znaleziony!")

        self.pdb_file = pdb_file

        # Wczytaj receptor za pomocą OpenBabel
        self.obconv = openbabel.OBConversion()
        self.obconv.SetInAndOutFormats("pdb", "pdb")

        self.mol = openbabel.OBMol()
        if not self.obconv.ReadFile(self.mol, pdb_file):
            raise RuntimeError("Nie udało się wczytać pliku PDB!")

    def add_hydrogens(self):
        """
        Dodaje brakujące atomy wodoru do receptora.
        """
        initial_atoms_count = self.mol.NumAtoms()
        self.mol.AddHydrogens()
        new_atoms_count = self.mol.NumAtoms()
        print(f"Dodano atomy wodoru. Liczba atomów wzrosła z {initial_atoms_count} do {new_atoms_count}.", flush=True)


    def minimize_with_constraints(self, stage1_iters=500, stage2_iters=50, tolerance=1.0e-4):
        """
        Przeprowadza minimalizację energii w dwóch etapach.

        Liczba iteracji dla etapu 1 (optymalizacja tylko wodorów przy więzach na ciężkich atomach)
        Liczba iteracji dla etapu 2 (łagodna minimalizacja całego układu)
        Tolerancja minimalizacji (jednostka domyślna zależna od pola siłowego)
        """
        print(f"Liczba atomów w molekule: {self.mol.NumAtoms()}", flush=True)

        # Użycie pola siłowego MMFF94
        ff = openbabel.OBForceField.FindForceField("MMFF94")
        if ff is None:
            raise RuntimeError("Nie znaleziono pola siłowego MMFF94!")

        if not ff.Setup(self.mol):
            raise RuntimeError("Nie udało się przygotować pola siłowego!")

        # Krok 1 – zamrożenie ciężkich atomów (atomic number > 1)
        for atom in openbabel.OBMolAtomIter(self.mol):
            if atom.GetAtomicNum() > 1:
                ff.AddFixedAtom(atom.GetIdx() - 1)  # Uwaga! OBForceField używa indeksów 0-based

        print("Rozpoczynam stage 1 minimalizacji (tylko wodory)...")
        ff.ConjugateGradients(stage1_iters, tolerance)
        ff.GetCoordinates(self.mol)

        # Krok 2 – łagodna minimalizacja całego układu
        ff.ClearFixedAtoms()
        print("Rozpoczynam stage 2 minimalizacji (cały układ)...")
        ff.ConjugateGradients(stage2_iters, tolerance)
        ff.GetCoordinates(self.mol)

        print("Minimalizacja zakończona.")

    def write_pdb(self, output_file):
        """
        Zapisuje zminimalizowany receptor do pliku PDB.
        """
        if not self.obconv.WriteFile(self.mol, output_file):
            raise RuntimeError(f"Nie udało się zapisać pliku {output_file}!")
        self.obconv.CloseOutFile()  # ważne, by zamknąć plik poprawnie
        print(f"Zminimalizowana struktura zapisana jako {output_file}")


if __name__ == '__main__':

    # Przykładowe ścieżki – dostosuj do swojej struktury
    input_pdb = "pdbqt_files/1hsg_receptor.pdb"              # Plik wejściowy PDB receptora (bez wodorów)
    output_pdb = "output_files/1hsg_receptor_minimized.pdb"  # Plik wynikowy z dodanymi wodorem i minimalizacją

    # Upewnij się, że folder output_files istnieje
    os.makedirs("output_files", exist_ok=True)

    # Utwórz obiekt minimalizatora, dodaj wodory i przeprowadź minimalizację
    minimizer = ReceptorMinimizerOpenBabel(input_pdb)
    minimizer.add_hydrogens()
    minimizer.minimize_with_constraints(stage1_iters=500, stage2_iters=50, tolerance=1.0e-4)
    minimizer.write_pdb(output_pdb)
