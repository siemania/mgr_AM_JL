#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Skrypt do minimalizacji receptora:
- Wczytuje receptor w formacie PDB (bez ligandów/wody)
- Dodaje atomy wodoru (przy pomocy RDKit)
- Minimalizuje receptor w dwóch krokach:
   1. Minimalizacja jedynie pozycji atomów wodoru (ciężkie atomy pozostają niezmienione)
   2. Krótka minimalizacja całego układu, by pozwolić na niewielkie przemieszczenia ciężkich atomów
- Zapisuje wynikowy receptor do nowego pliku PDB
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os

class ReceptorMinimizer:
    def __init__(self, pdb_file):
        """
        Inicjalizacja – wczytanie receptora z pliku PDB.
        """
        self.added_h_count = None
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"Plik {pdb_file} nie został znaleziony!")
        # Wczytujemy receptor – zakładamy, że w PDB nie ma hydrogenu (czysty receptor)
        self.mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if self.mol is None:
            raise ValueError(f"Nie udało się wczytać pliku PDB: {pdb_file}")

    def add_hydrogens(self):
        """
        Dodaje wszystkie niezbędne atomy wodoru do receptora oraz liczy ich liczbę.
        """
        initial_h_count = sum(1 for atom in self.mol.GetAtoms() if atom.GetAtomicNum() == 1)
        self.mol = Chem.AddHs(self.mol)
        final_h_count = sum(1 for atom in self.mol.GetAtoms() if atom.GetAtomicNum() == 1)
        self.added_h_count = final_h_count - initial_h_count

    def minimize_with_constraints(self, free_iters=50, mmff_max_iters=200):
        """
        Przeprowadza dwustopniową minimalizację:
          - Krok 1: Minimalizacja z zablokowanymi ciężkimi atomami,
                    czyli optymalizacja pozycji wyłącznie dodanych atomów wodoru.
                    Używamy tu standardowego optymalizatora MMFF, przekazując listę indeksów ciężkich atomów
          - Krok 2: Minimalizacja całego układu przez ograniczoną liczbę iteracji (free_iters),
                    aby pozwolić na niewielkie przesunięcia ciężkich atomów.
        """
        # Wybieramy indeksy ciężkich atomów – atomy o liczbie atomowej większej niż 1
        heavy_atom_indices = [atom.GetIdx() for atom in self.mol.GetAtoms() if atom.GetAtomicNum() > 1]

        # Krok 1: Minimalizacja H (ciężkie atomy zablokowane)
        if heavy_atom_indices:
            # Dla MMFF, argument fixedAtoms zamraża podane atomy, optymalizując pozostałe.
            # W ten sposób pozycje ciężkich atomów pozostaną bez zmian.
            success1 = AllChem.MMFFOptimizeMolecule(self.mol, maxIters=mmff_max_iters, fixedAtoms=heavy_atom_indices)
            if success1 != 0:
                print("Ostrzeżenie: minimalizacja atomów wodoru mogła nie zakończyć się optymalnie.")

        # Krok 2: Krótka minimalizacja całego układu – lekkie przesunięcia ciężkich atomów
        success2 = AllChem.MMFFOptimizeMolecule(self.mol, maxIters=free_iters)
        if success2 != 0:
            print("Ostrzeżenie: końcowa minimalizacja mogła nie zakończyć się optymalnie.")

    def write_pdb(self, output_file):
        """
        Zapisuje zminimalizowany receptor do pliku PDB.
        """
        writer = Chem.PDBWriter(output_file)
        writer.write(self.mol)
        writer.close()


if __name__ == '__main__':
    # Przykładowe ścieżki – dostosuj do swojej struktury katalogów
    input_pdb = "pdb_files\\1hsg_receptor.pdb"            # Ścieżka do wejściowego receptor.pdb
    output_pdb = "output_files\\receptor_1hsg_minimized.pdb"  # Ścieżka do pliku wynikowego

    # Upewnij się, że katalog output_files istnieje
    os.makedirs("output_files", exist_ok=True)

    # Utwórz obiekt minimalizacji
    minimizer_file = ReceptorMinimizer(input_pdb)

    # Dodaj atomy wodoru
    print("Dodawanie atomów wodoru...", flush=True)
    minimizer_file.add_hydrogens()
    print(f"Atomów wodoru dodano: {minimizer_file.added_h_count}", flush=True)

    # Wykonaj minimalizację z restrykcjami
    print("Rozpoczynam minimalizację receptora...", flush=True)
    # minimizer_file.minimize_with_constraints(free_iters=50)
    print("Minimalizacja zakończona.", flush=True)

    # Zapisz wynik do pliku
    minimizer_file.write_pdb(output_pdb)
    print(f"Zminimalizowany receptor zapisany w: {output_pdb}", flush=True)
