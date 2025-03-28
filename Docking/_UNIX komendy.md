
### W LINUX wyłapuje 1 od góry "HET" z PDB i oznacza jako liganda, następnie upewnia się w PDB że wybieramy 1 model (jeśli jest ich kilka) i wyszukuje po hetID liganda zapisując go do pliku
hetID=$(awk '/\<^HET\>/ {print $2}' 1ttv.pdb | awk 'NR==1'); awk '/MODEL        1/{flag=1} flag; /ENDMDL/{flag=0}' 1ttv.pdb | grep "$hetID" | grep '^HETATM' > 1ttv_ligand.pdb

# przypadek kiedy nie ma wielu modeli
hetID=$(awk '/\<^HET\>/ {print $2}' 1hsg.pdb | awk 'NR==1'); grep "$hetID" 1hsg.pdb | grep '^HETATM'