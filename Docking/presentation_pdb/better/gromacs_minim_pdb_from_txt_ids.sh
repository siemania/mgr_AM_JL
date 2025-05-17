#!/bin/bash

while read -r id _; do
    gmx pdb2gmx -f ${id}_reduced.pdb \
        -o ${id}_reduced.gro \
        -p topol.top \
        -i posre.itp \
        -ff oplsaa \
        -water none \
        -ignh && \
    gmx editconf -f ${id}_reduced.gro -o ${id}_box.gro -c -d 1.0 && \
    gmx grompp -f minim.mdp -c ${id}_box.gro -p topol.top -maxwarn 1 -o ${id}_minim.tpr && \
    gmx mdrun -v -deffnm ${id}_minim && \
    echo 0 | gmx trjconv -f ${id}_minim.gro -s ${id}_minim.tpr -pbc mol -ur compact -o ${id}_no_pbc.gro && \
    gmx editconf -f ${id}_no_pbc.gro -o ${id}_minimized.pdb
done < ../presentation_pdb.txt

