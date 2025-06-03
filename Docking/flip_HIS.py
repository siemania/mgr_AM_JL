import pymol2
from math import radians

def flip_HIS(pdb_path, output_path, residue_id):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb_path, 'prot')

        # Definiowanie atomow torsyjnych (tworza kat dwuscienny)
        sele = f"(resi {residue_id} and name CA+CB+CG+ND1)"
        atoms = pymol.cmd.get_model(sele).atom

        if len(atoms) == 4:
            atom_names = [a.name for a in atoms]
            ids = [a.index for a in atoms]

            pymol.cmd.set_dihedral(
                f"id {ids[0]}",
                f"id {ids[1]}",
                f"id {ids[2]}",
                f"id {ids[3]}",
                180 # obrot o 180 stopni
            )
            print(f"Wykonano flip HIS:{residue_id}.")
        else:
            print(f"Błąd: Nie znaleziono pełnego zestawu atomów w reszcie {residue_id}.")

        pymol.cmd.save(output_path, 'prot')