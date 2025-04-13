import subprocess


def minimize_receptor(input_pdb, output_pdb):
    """
    Minimizes the receptor using Open Babel. The input is a PDB file, and the output
    is the minimized PDB structure.

    :param input_pdb: The PDB file of the receptor
    :param output_pdb: The output minimized PDB file
    """
    # Perform energy minimization with Open Babel
    try:
        subprocess.run(["obabel", input_pdb, "-O", output_pdb, "--minimize", "--ff", "ghemical"], check=True)
        print(f"Receptor minimized and saved as {output_pdb}")
    except subprocess.CalledProcessError:
        print("Error: Minimization failed.")
        return None

    return output_pdb


def prepare_receptor_for_autodock(minimized_pdb, output_pdbqt):
    """
    Prepare the minimized receptor file for AutoDock (convert to PDBQT format)

    :param minimized_pdb: The minimized PDB file
    :param output_pdbqt: The output PDBQT file for AutoDock
    """
    try:
        # Use AutoDockTools (ADT) to prepare the receptor for docking
        subprocess.run(["prepare_receptor4.py", "-r", minimized_pdb, "-o", output_pdbqt], check=True)
        print(f"Receptor prepared for AutoDock and saved as {output_pdbqt}")
    except subprocess.CalledProcessError:
        print("Error: Unable to prepare the receptor file for AutoDock.")
        return None

    return output_pdbqt


# Main execution block
if __name__ == "__main__":
    # Input and output files
    input_receptor_pdb = "receptor.pdb"  # Input receptor PDB file
    minimized_receptor_pdb = "receptor_minimized.pdb"  # Output minimized receptor PDB
    receptor_pdbqt = "receptor_minimized.pdbqt"  # Receptor in PDBQT format for AutoDock

    # Step 1: Minimize the receptor structure
    minimized_receptor = minimize_receptor(input_receptor_pdb, minimized_receptor_pdb)

    if minimized_receptor:
        # Step 2: Prepare the receptor for AutoDock
        prepare_receptor_for_autodock(minimized_receptor, receptor_pdbqt)
