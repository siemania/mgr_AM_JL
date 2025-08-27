

def modify_gpf(input_gpf, output_gpf, map_dir="map_grid_files", receptor_dir="pdbqt_files", quiet=False):
    """Modyfikuje i zapisuje w nowej lokalizacji plik .gpf uwzględniając inne foldery dla map energetycznych
    oraz receptora"""
    with open(input_gpf, "r") as f:
        lines = f.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith("map "):
            # Pobiera nazwę pliku mapy: map 1hsg_receptor.A.map upewniając się że pobiera tylko tę nazwę
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"map {map_dir}/{filename}              # atom-specific affinity map\n")
        elif line.startswith("elecmap "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"elecmap {map_dir}/{filename}          # electrostatic potential map\n")
        elif line.startswith("dsolvmap "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"dsolvmap {map_dir}/{filename}        # desolvation potential map\n")
        elif line.startswith("receptor "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"receptor {receptor_dir}/{filename}         # macromolecule\n")
        elif line.startswith("gridfld "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"gridfld {map_dir}/{filename}       # grid_data_file\n")
        elif line.startswith("fld "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"fld {map_dir}/{filename}           # grid_data_file\n")
        elif line.startswith("desolvmap "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"desolvmap {map_dir}/{filename}        # desolvation map\n")
        elif line.startswith("move "):
            filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
            modified_lines.append(f"move {receptor_dir}/{filename}               # small molecule\n")
        else:
            modified_lines.append(line) # Zostawiamy linię bez zmian

    with open(output_gpf, "w") as f:
        f.writelines(modified_lines)

    if quiet != True: print(f"Plik {output_gpf} został zapisany.")
    return


def modify_gdpf_overwrite(input_gpf, map_dir="map_grid_files", receptor_dir="pdbqt_files", quiet=False):
    """Modyfikuje i nadpisuje plik .gpf uwzględniając inne foldery dla map energetycznych oraz receptora"""
    # Wczytanie pliku i sprawdzenie, czy już zawiera poprawne ścieżki
    with open(input_gpf, "r+") as f:
        lines = f.readlines()
        f.seek(0)  # Cofnięcie wskaźnika na początek pliku
        f.truncate()  # Usunięcie starej zawartości

        for line in lines:
            if line.startswith("map "):
                # Pobiera nazwę pliku mapy: map 1hsg_receptor.A.map upewniając się że pobiera tylko tę nazwę
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"map {map_dir}/{filename}              # atom-specific affinity map\n"
            elif line.startswith("elecmap "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"elecmap {map_dir}/{filename}          # electrostatic potential map\n"
            elif line.startswith("dsolvmap "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"dsolvmap {map_dir}/{filename}        # desolvation potential map\n"
            elif line.startswith("receptor "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"receptor {receptor_dir}/{filename}         # macromolecule\n"
            elif line.startswith("gridfld "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"gridfld {map_dir}/{filename}       # grid_data_file\n"
            elif line.startswith("fld "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"fld {map_dir}/{filename}           # grid_data_file\n"
            elif line.startswith("desolvmap "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"desolvmap {map_dir}/{filename}        # desolvation map\n"
            elif line.startswith("move "):
                filename = line.split()[1].rsplit(sep="/", maxsplit=1)[-1]
                new_line = f"move {receptor_dir}/{filename}               # small molecule\n"
            else:
                new_line = line  # Zostawiamy linię bez zmian

            # Zapisujemy przetworzoną linię do pliku
            f.write(new_line)

    if quiet != True: print(f"Plik {input_gpf} został zaktualizowany.") 
    return


def modify_pdbqt_overwrite(input_pdbqt, quiet=False):
    """Modyfikuje i nadpisuje plik .pdbqt usuwając wiersze na podstawie alternatywnej pozycji z kolumny 17 w pliku
    receptora, gdzie znajduje się znak inny niż A lub biały znak, końcowo pozostałe wiersze w kolumnie 17 zastąpi znak
    biały"""
    # Wczytanie pliku i sprawdzenie, czy już zawiera poprawne ścieżki
    with open(input_pdbqt, "r+") as f:
        lines = f.readlines()
        f.seek(0)  # Cofnięcie wskaźnika na początek pliku
        f.truncate()  # Usunięcie starej zawartości

        modified_lines = []
        for line in lines:
            if line[20] != " ": # Jeżeli znak 21 w PDB zawiera jakiś znak
                temp = line[17] # Zachowujemy znak z kolumny 18 w PDB
                new_line = line[:17] + "" + line[18:] # Usuwamy znak z kolumny 18 w PDB
                new_line = new_line[:12] + temp + new_line[13:] # Dodajemy zachowany znak do kolumny 13 tak jak robi to ADT
            else:
                new_line = line  # Zostawiamy linię bez zmian

            # Zapisujemy przetworzoną linię do pliku
            f.write(new_line)

        if quiet != True: print(f"Plik {input_pdbqt} został zaktualizowany.")
        return


def modify_fld_overwrite(input_fld, quiet=False):
    """Modyfikuje i nadpisuje plik .fld usuwając wystąpienia 'map_grid_files' ze ścieżek"""
    with open(input_fld, "r+") as f:
        content = f.read()
        f.seek(0)
        f.truncate()
        modified_content = content.replace("map_grid_files/", "")
        f.write(modified_content)

    if quiet != True: print(f"Plik {input_fld} został zaktualizowany.")
    return

if __name__ == '__main__':
    # Przykładowe użycie
    # modify_gdpf("grid_dock_files/1hsg_grid.gpf", "output.gpf")
    modify_gdpf_overwrite("grid_dock_files/1hsg_grid.gpf")
    # modify_fld_overwrite("map_grid_files/example.fld")
