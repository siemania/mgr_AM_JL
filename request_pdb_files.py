import os
import requests
import pandas as pd

def download_pdb_files(csv_file, folder_pdb, folder_ent):
    """Pobiera pliki PDB na podstawie ID z pliku CSV."""
	
    # Wczytanie pliku CSV
    df = pd.read_csv(csv_file)
    
    # Sprawdzenie czy kolumna "ID" istnieje
    if "id" not in df.columns:
        raise ValueError("Brak kolumny 'ID' w pliku CSV")
    
    # Tworzenie folderow na pliki
    os.makedirs(folder_pdb, exist_ok=True)
    # os.makedirs(folder_ent, exist_ok=True)

    base_url = "https://files.rcsb.org/download/{}.pdb"
    missing_files = []
    
    for pdb_id in df["id"].dropna().astype(str):
        pdb_id_clean = pdb_id.strip().lower()  # Usunięcie białych znaków


        pdb_path = os.path.join(folder_pdb, f"{pdb_id_clean}.pdb")
        # ent_path = os.path.join(folder_ent, f"{pdb_id_clean}.entOrg")

        # Pobieranie plikow pdb
        try:
            response = requests.get(base_url.format(pdb_id_clean), timeout=10)
            if response.status_code == 200: # Działa
                with open(pdb_path, "wb") as f:
                    f.write(response.content)
                print(f"Pobrano: {pdb_id_clean}.pdb")

               # Skopiowanie pliku z rozszerzeniem .ent
               # with open (ent_path, "wb") as f:
               #    f.write(response.content)
               # print(f"Zapisano jako: pdb{pdb_id_clean}.entOrg")
            else:
                print(f"Nie udało się pobrać: {pdb_id_clean}")
                missing_files.append(pdb_id_clean)
                
        except requests.RequestException as e:
            print(f"Błąd przy pobieraniu {pdb_id_clean}.pdb: {e}")
            missing_files.append(pdb_id_clean)
    
    # Podsumowanie
    if missing_files:
        print("\nNie udało się pobrać następujących plików:")
        for m in missing_files:
            print(" -", m)
    else:
        print("\nWszystkie pliki zostały pobrane i zapisane pomyślnie.")



if __name__ == "__main__":
    
    csv_file = "baza_ids.csv"  # Plik CSV z ID struktur
    folder_pdb = "pdb_files"  # Folder na pobrane pliki
    folder_ent = "ent_files"
    
    download_pdb_files(csv_file, folder_pdb, folder_ent)