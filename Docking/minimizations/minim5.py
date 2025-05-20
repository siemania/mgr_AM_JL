from rdkit import Chem
from rdkit.Chem import AllChem

# Wczytaj plik PDB ze współrzędnymi
mol = Chem.MolFromPDBFile("pdb1crn.ent", removeHs=False)

if mol is None:
    raise ValueError("Nie udało się wczytać pliku PDB.")

# Jeśli są brakujące wodory – dodaj
mol = Chem.AddHs(mol, addCoords=True)  # dodaje H z wyznaczeniem pozycji

# Nie ruszamy geometrii! (brak EmbedMolecule)

# Zbuduj pole sił
props = AllChem.MMFFGetMoleculeProperties(mol)
ff = AllChem.MMFFGetMoleculeForceField(mol, props)

# Zamroź atomy ciężkie
for i, atom in enumerate(mol.GetAtoms()):
    if atom.GetAtomicNum() != 1:  # nie-H
        ff.AddFixedPoint(i)

# Minimalizuj tylko wodory
ff.Minimize()

# Zapisz wynik
writer = Chem.rdmolfiles.PDBWriter("ligand_minimized_preserved_geometry.pdb")
writer.write(mol)
writer.close()
