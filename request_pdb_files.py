import os
import requests
import pandas as pd

def download_pdb_files(csv_file, output_folder):
    """Pobiera pliki PDB na podstawie ID z pliku CSV."""
    # Wczytanie pliku CSV
    df = pd.read_csv(csv_file)
    
    # Sprawdzenie czy kolumna "ID" istnieje
    if "id" not in df.columns:
        raise ValueError("Brak kolumny 'ID' w pliku CSV")
    
    # Tworzenie folderu na pliki
    os.makedirs(output_folder, exist_ok=True)
    
    base_url = "https://files.rcsb.org/download/{}.pdb"
    missing_files = []
    
    for pdb_id in df["id"].dropna().astype(str):
        pdb_id = pdb_id.strip()  # Usunięcie białych znaków
        url = base_url.format(pdb_id)
        file_path = os.path.join(output_folder, f"{pdb_id}.pdb")
        
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200: # Działa
                with open(file_path, "wb") as f:
                    f.write(response.content)
                print(f"Pobrano: {pdb_id}")
            else:
                print(f"Nie udało się pobrać: {pdb_id}")
                missing_files.append(pdb_id)
                
        except requests.RequestException as e:
            print(f"Błąd przy pobieraniu {pdb_id}: {e}")
            missing_files.append(pdb_id)
    
    # Podsumowanie
    if missing_files:
        print("Nie udało się pobrać następujących plików:")
        for pdb in missing_files:
            print(pdb)
    else:
        print("Wszystkie pliki zostały pobrane pomyślnie.")



if __name__ == "__main__":
    
    csv_file = "baza_ids.csv"  # Plik CSV z ID struktur
    output_folder = "pdb_files"  # Folder na pobrane pliki
    
    download_pdb_files(csv_file, output_folder)