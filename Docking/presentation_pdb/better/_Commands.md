# pdbfixer
for f in ../*.pdb; do pdbfixer "${f}" --verbose --add-atoms=all --keep-heterogens=none --replace-nonstandard --add-residues --output="${f%.pdb}_fixed.pdb"; done

# reduce
for f in *_fixed.pdb; do reduce -BUILD ${f} > ${f%_fixed.pdb}_reduced.pdb; done

# Variable
file="ID TWOJEGO PDB (BEZ CUDZYSŁOWIU)"

# GROMACS
gmx pdb2gmx -f ${file}_reduced.pdb \
    -o ${file}_reduced.gro \
    -p topol.top \
    -i posre.itp \
    -ff oplsaa \
    -water none \
    -ignh

# Tworzenie boxa
gmx editconf -f ${file}_reduced.gro -o ${file}_box.gro -c -d 1.0

# Plik minim.mdp
integrator  = steep
nsteps      = 50000
emtol       = 1000.0
emstep      = 0.01
nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz

# Generowanie pliku .tpr
gmx grompp -f minim.mdp -c ${file}_box.gro -p topol.top -maxwarn 1 -o ${file}_minim.tpr

# Uruchomienie minimalizacji
gmx mdrun -v -deffnm ${file}_minim

# Przywraca 
gmx trjconv -f ${id}_minim.gro -pbc mol -ur compact -o ${id}_no_pbc.gro

# Konwersja do PDB
gmx editconf -f ${file}_minim.gro -o ${file}_minimized.pdb

# Przywrócenie jak najbliżej systemu koordynatów
echo 0 | gmx trjconv -f twoja_struktura.pdb -s oryginalna_referencja.pdb -fit rot+trans -o dopasowana_struktura.pdb

# Przygotowanie pliku PDBQT (w python i reszta) [Robiłem ręcznie w AutoDock, a potem automatico]
# Merge non-polar, calculate Kollman charges, assign AD4 type, grid > choose macromolecule > save .pdbqt
prepare_receptor4.py -r ${file}_minimized.pdb \
    -o ${file}_receptor.pdbqt \
    -U nphs_lps_waters \

