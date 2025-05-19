
# W LINUX wyłapuje 1 od góry "HET" z PDB i oznacza jako liganda, następnie upewnia się w PDB że wybieramy 1 model (jeśli jest ich kilka) i wyszukuje po hetID liganda zapisując go do pliku
hetID=$(awk '/\<^HET\>/ {print $2}' 1ttv.pdb | awk 'NR==1'); awk '/MODEL        1/{flag=1} flag; /ENDMDL/{flag=0}' 1ttv.pdb | grep "$hetID" | grep '^HETATM' > 1ttv_ligand.pdb

# przypadek kiedy nie ma wielu modeli
hetID=$(awk '/\<^HET\>/ {print $2}' 1hsg.pdb | awk 'NR==1'); grep "$hetID" 1hsg.pdb | grep '^HETATM'

# Wyświetla ostatnie linijki co sekundę końcowego dokowania, tylko jeśli się nie zakończyły
watch -n 1 'for f in *.dlg; do tail -n 3 "$f" | grep -q "^___" && continue; echo "==> $f <=="; tail -n 3 "$f"; done'

# Modyfikacja parametr extension
${variable#pattern} – usuwa najkrótszy fragment pasujący do wzorca z początku wartości.
${variable##pattern} – usuwa najdłuższy fragment pasujący do wzorca z początku.
${variable%pattern} – usuwa najkrótszy fragment pasujący do wzorca z końca.
${variable%%pattern} – usuwa najdłuższy fragment pasujący do wzorca z końca.
${variable/default} – jeśli zmienna jest pusta, podstawia wartość domyślną.

# Przykład dla .pdbq to .pdb
for file in *.pdbq; do mv "$file" "${file%.pdbq}.pdb"; done

# Kopiowanie plików na podstawie pliku tekstowego
while read -r id _; do cp "pdb_files/${id}.pdb" presentation_pdb/; done < presentation_pdb.txt