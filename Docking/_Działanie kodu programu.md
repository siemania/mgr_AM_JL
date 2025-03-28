# =============================================================================
#                       Właściwy program
# =============================================================================
# Przygotowanie receptora bezpośrednio z PDB (PDB → PDBQT)
# Dodawane są ładunki Kollmana
# -A "bonds_hydrogens" → build a single bond from each atom with no bonds to its closest neighbor and add hydrogens
# -U "nphs" → merge charges and remove non-polar hydrogens; "waters": remove water residues; "lps": merge charges and remove lone pairs
# -e "True" → delete every nonstd residue from any chain
run_command(f'"{python}" prepare_receptor4.py -r {receptor_pdb} -o {receptor_pdbqt} -A "bonds_hydrogens" -e "True" -U "nphs" "waters" "lps"')
print("Receptor gotowy")

# Przygotowanie liganda (PDB → PDBQT)
# Dodawane są ładunki Gasteigera
run_command(f'"{python}" prepare_ligand4.py -l {ligand_pdb} -o {ligand_pdbqt}')
print("Ligand gotowy")

# Przygotowanie pliku parametrów siatki (.GPF)
#!!! Na ten moment wszystkie parametry są takie jakie są
run_command(f'"{python}" prepare_gpf4.py -l {ligand_pdbqt} -r {receptor_pdbqt} -o {gpf}')
modify_gdpf_overwrite(f"grid_dock_files/{file_name}_grid.gpf", quiet=True) # Poprawki lokalizacyjne w pliku
print("Grid gotowy")

# Uruchomienie AutoGrid (GPF → GLG plik typu log)
run_command(f'{autogrid} -p {gpf} -l {glg}')
print("Mapy energetyczne gotowe")

# Przygotowanie pliku parametrów dokowania (.DPF)
#!!! Na ten moment wszystkie parametry są takie jakie są
run_command(f'"{python}" prepare_dpf42.py -l {ligand_pdbqt} -r {receptor_pdbqt} -o {dpf}')
modify_gdpf_overwrite(f"grid_dock_files/{file_name}_dock.dpf", quiet=True) # Poprawki lokalizacyjne w pliku
print("Parametry do dokowania gotowe")

# Uruchomienie AutoDock (DPF → DLG plik typu log)
run_command(f'{autodock} -p {dpf} -l {dlg}')
print("Dokowanie zakończone pomyślnie")

# Zapisanie najlepszego kompleksu
run_command(f'"{python}" write_all_complexes.py -d {dlg} -r {receptor_pdbqt} -o {file_name}_bestcomplex -b')
print("Utworzono najlepszy kompleks ligand-receptor")