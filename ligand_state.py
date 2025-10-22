from os.path import splitext
import argparse

def main(args):
    
    filename = splitext(*args.file)[0]
    
    with open(f"{filename}.dlg", encoding="utf-8") as file:
        
        if args.number:
            run_docked = int(args.number)
        else:
            # Find best docked state
            for line in file:
                if line.startswith("Rank | Energy"): # End of clustering
                    next(file)
                    run_docked = int(next(file)[18:23])
                    break
                
        file.seek(0)
        
        occ_temp = []
        for line in file:

            # Check for best ligand docked    
            if line.startswith("DOCKED: USER    Run"):
                run = int(line.split()[-1])
                if run == run_docked:
                    break                    
        
            # Check for input ligand to copy 2 columns
            elif line.startswith("INPUT-LIGAND-PDBQT: HETATM"):
                occ_temp.append(line[74:86])

        # Continue getting ligand
        ligand = []
        charges = []
        for line in file:
            if line.startswith(("DOCKED: HETATM", "DOCKED: ATOM")): # AD4
                ligand.append(line[8:62])
                charges.append(line[74:88])
                
    # Concat lists
    union = []
    for part1, part2, part3 in zip(ligand, occ_temp, charges):
        u = ''.join((part1, part2, part3, "\n"))
        union.append(u)
    
    # Save ligand
    with open(f"{filename}.pdb", "w", encoding="utf-8") as f:
        f.writelines(union)
                


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog="Ligand getter from docking states",
        description="From dlg files it gets ligand based on chosen run number",
        usage="Usage: python ligand_state.py -f filename.dlg [opt: -n number]"
    )
    
    parser.add_argument("-f", "--file", nargs=1, help="Input .dlg file")
    parser.add_argument("-n", "--number", nargs="?", help="Run number for ligand",
                        default=None)
    
    args = parser.parse_args()
    
    main(args)
