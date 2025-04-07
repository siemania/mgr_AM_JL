# -*- coding: utf-8 -*-
import os
import sys
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem



# Scieżka do folderu z plikami receptorów z dodanymi wodorami
username = os.environ.get('USERNAME') or os.environ.get('USER')
PROJECT_ROOT = "/home/" + username + "/IdeaProjects/mgr_AM_JL"
TEST_SKRYPTU_DIR = os.path.join(PROJECT_ROOT, "Docking/TEST_SKRYPTU")
PDB_FILE_DIR = os.path.join(TEST_SKRYPTU_DIR, "pdb_files")
MINIMIZATION_FILE_DIR = os.path.join(TEST_SKRYPTU_DIR, "minimization")
MINI1_FILE_DIR = os.path.join(MINIMIZATION_FILE_DIR, "mini1")
MINI2_FILE_DIR = os.path.join(MINIMIZATION_FILE_DIR, "mini2")

def huge_minimization(pdb_file):
    """Ogolna minimalizacja wszystkiego."""
    # Zaladowanie pliku
    print(pdb_file)
    file_converter = openbabel.OBConversion()
    file_converter.SetInAndOutFormats("pdb", "pdb")
    # Tworzenie nowej molekuly
    mol = openbabel.OBMol()
    success = file_converter.ReadFile(mol, pdb_file)
    if not success:
        raise ValueError ("Nie udalo się wczytac pliku PDB.")


    mol.SetAutomaticFormalCharge(False)
    mol.ConnectTheDots()
    mol.FindRingAtomsAndBonds()
    mol.PerceiveBondOrders()
    #mol.ClearAromaticFlags()
    mol.AddHydrogens()

    # Ustawianie pola silowego
    force_field = openbabel.OBForceField.FindForceField("UFF")
    if force_field is None:
        raise ValueError("Nie znaleziono pola silowego UFF.")
    force_field.Setup(mol)
    force_field.SteepestDescent(10)
    force_field.GetCoordinates(mol)

    # Zapisanie zminimalizowanej struktury
    output_file = os.path.basename(pdb_file).split(".")[0]
    output_file = os.path.join(MINI1_FILE_DIR, output_file)
    file_converter.WriteFile(mol, output_file)

    return output_file


def minimization_with_constrains (output_file_1):
    """Mniejsza minimalizacja z nalozeniem wiezow."""
    # Zaladowanie pliku PDB po pierwszej minimalizacji
    mol = Chem.MolFromPDBFile(output_file_1, removeHs = False)

    # Sprawdzenie czy struktura zostala poprawnie wczytana
    if mol is None:
        raise ValueError ("Nie udalo się wczytac pliku PDB.")

    # Ustawienie pola silowego
    force_field = AllChem.UFFGetMoleculeForceField(mol)

    # Tworzenie listy z indeksami atomow ciezkich i zakladanie wiezow
    heavy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    for index in heavy_atoms:
        force_field.AddFixedPoint(index)

    # Przeprowadzenie minimalizacji struktury z zalozonymi wiezami i zapisanie zminimalizowanej struktury
    force_field.Minimize()
    output_file_2 = os.path.basename(output_file_1).split(".")[0]
    output_file_2 = os.path.join(MINI2_FILE_DIR, output_file_2)
    Chem.MolToPDBFile(mol, output_file_2)

    return output_file_2

def optimize_rec (PDB_FILE_DIR):
    """Obsluga argumentow wywolania skryptu."""
    for pdb_file in os.listdir(PDB_FILE_DIR):
        print(PDB_FILE_DIR, pdb_file)
        if not pdb_file.endswith(".pdb"):
            continue
        else:
            pdb_path = os.path.join(PDB_FILE_DIR, pdb_file)
            output_file_1 = huge_minimization(pdb_path)
            minimization_with_constrains(output_file_1)



if __name__ == "__main__":
    #if len(sys.argv) < 2:
    #    print ("Brak sciezki do folderu")
    #    sys.exit(1)
    #else:

    #PDB_FILE_DIR = sys.argv[1]
    for path in [MINIMIZATION_FILE_DIR, MINI1_FILE_DIR, MINI2_FILE_DIR]:
        if not os.path.isdir(path):
            os.makedirs(path)

    optimize_rec(PDB_FILE_DIR)

