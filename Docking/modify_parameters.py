# -*- coding: utf-8 -*-
import os

def modify_gpf(input_gpf, output_gpf, map_dir="map_grid_files", receptor_dir="pdbqt_files", quiet=False):
    """Modyfikuje i zapisuje w nowej lokalizacji plik .gpf uwzględniając inne foldery dla map energetycznych
    oraz receptora"""
    with open(input_gpf, "r") as f:
        lines = f.readlines()

    modified_lines = []
    for line in lines:
        if line.startswith("map "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("map %s/%s              # atom-specific affinity map\n" % (map_dir, filename))
        elif line.startswith("elecmap "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("elecmap %s/%s          # electrostatic potential map\n" % (map_dir, filename))
        elif line.startswith("dsolvmap "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("dsolvmap %s/%s        # desolvation potential map\n" % (map_dir, filename))
        elif line.startswith("receptor "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("receptor %s/%s         # macromolecule\n" % (receptor_dir, filename))
        elif line.startswith("gridfld "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("gridfld %s/%s       # grid_data_file\n" % (map_dir, filename))
        elif line.startswith("fld "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("fld %s/%s           # grid_data_file\n" % (map_dir, filename))
        elif line.startswith("desolvmap "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("desolvmap %s/%s        # desolvation map\n" % (map_dir, filename))
        elif line.startswith("move "):
            filename = line.split()[1].rsplit("/", 1)[-1]
            modified_lines.append("move %s/%s               # small molecule\n" % (receptor_dir, filename))
        else:
            modified_lines.append(line)

    with open(output_gpf, "w") as f:
        f.writelines(modified_lines)

    if not quiet:
        print("Plik %s został zapisany." % output_gpf)
    return


def modify_gdpf_overwrite(input_gpf, map_dir="map_grid_files", receptor_dir="pdbqt_files", quiet=False):
    """Modyfikuje i nadpisuje plik .gpf uwzględniając inne foldery dla map energetycznych oraz receptora"""
    with open(input_gpf, "r+") as f:
        lines = f.readlines()
        f.seek(0)
        f.truncate()

        for line in lines:
            if line.startswith("map "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "map %s/%s              # atom-specific affinity map\n" % (map_dir, filename)
            elif line.startswith("elecmap "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "elecmap %s/%s          # electrostatic potential map\n" % (map_dir, filename)
            elif line.startswith("dsolvmap "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "dsolvmap %s/%s        # desolvation potential map\n" % (map_dir, filename)
            elif line.startswith("receptor "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "receptor %s/%s         # macromolecule\n" % (receptor_dir, filename)
            elif line.startswith("gridfld "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "gridfld %s/%s       # grid_data_file\n" % (map_dir, filename)
            elif line.startswith("fld "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "fld %s/%s           # grid_data_file\n" % (map_dir, filename)
            elif line.startswith("desolvmap "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "desolvmap %s/%s        # desolvation map\n" % (map_dir, filename)
            elif line.startswith("move "):
                filename = line.split()[1].rsplit("/", 1)[-1]
                new_line = "move %s/%s               # small molecule\n" % (receptor_dir, filename)
            else:
                new_line = line

            f.write(new_line)

    if not quiet:
        print("Plik %s został zaktualizowany." % input_gpf)
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


if __name__ == '__main__':
    modify_gdpf_overwrite("grid_dock_files/1hsg_grid.gpf")
