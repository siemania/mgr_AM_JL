import os

def extract_ligand(pdb_file, output_folder):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    models = False
    het_id = None
    max_atom_number = -1
    
    # Znalezienie HET z najwyższą wartością numeracji atomów (kolumny 21-25)
    for line in lines:
        if line.startswith("HET "):
            atom_number = int(line[20:25].strip())  # Pobranie numeru atomu (21-25 kolumna)
            current_het_id = line[7:10].strip()  # Pobranie het_id (8-10 kolumna)
            if atom_number > max_atom_number:
                max_atom_number = atom_number
                het_id = current_het_id
        if line.startswith("NUMMDL "):  # Sprawdzenie czy są modele
            models = True 
        if line.startswith("ATOM "):  # Pobieramy tylko pierwszy HET
            break
    
    if not het_id:
        print(f"Brak ligandu w {pdb_file}")
        return
    
    extracted_lines = []
    inside_model = False # Na początku brak modelu
    
    if models:        
        # Wyszukuje wszystkie wpisy w obrębie modelu 1
        for line in lines:
            if line.startswith("MODEL        1"):
                inside_model = True
            elif line.startswith("ENDMDL"):
                inside_model = False
                break
            if inside_model and line.startswith("HETATM") and het_id in line:
                extracted_lines.append(line)
    else:
        # Wyszukuje wpisy BEZ określonych modeli w PDB
        for line in lines:
            if line.startswith("HETATM") and het_id in line:
                extracted_lines.append(line)
    
    if extracted_lines:
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, f"{os.path.basename(pdb_file).replace('.pdb', '')}_ligand.pdb")
        log_line = f"Zapisano ligand {het_id} z {os.path.basename(pdb_file)} do: {output_file}"
        
        with open(output_file, 'w') as f:
            f.writelines(extracted_lines) # Zapis do pliku
        print(log_line)
        
    else:
        print(f"Brak odpowiednich wpisów HETATM w {pdb_file}")


if __name__ == '__main__':
    
    # Foldery wejściowy i wyjściowy
    input_folder = "pdb_files"
    output_folder = "ligands"
    
    # Przetwarzanie każdego pliku PDB
    for PDB_file in os.listdir(input_folder):
        if PDB_file.endswith(".pdb"):
            extract_ligand(os.path.join(input_folder, PDB_file), output_folder)
    input("~~~~~~~~Naciśnij dowolny przycisk by zakończyć~~~~~~~~")
