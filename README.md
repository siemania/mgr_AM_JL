# Automatyzacja Dokowania

Projekt **mgr_AM_JL** to zestaw skryptów mających na celu automatyzację przygotowania receptorów i ligandów do dokowania. Początkowo projekt automatycznie wykonuje serię poleceń związanych z przygotowaniem receptorów i ligandów przy użyciu MGLTools oraz AutoDock/AutoDock GPU, a w dalszej perspektywie planowane jest rozwinięcie narzędzi do:

- Optymalizacji konformacji wybranych reszt (np. HIS, GLN, ASN)
- Dokowania z wodą i bez wody
- Optymalizacji układu wodorów (obrót grup metylowych itp.)
- Dodatkowych algorytmów optymalizacyjnych wykorzystujących pola siłowe z programów jak np. OpenMM, OpenBabel, GROMACS

## Spis treści

- [Wymagania](#wymagania)
- [Instalacja](#instalacja)
- [Uruchomienie programu](#uruchomienie-programu)
- [Struktura projektu](#struktura-projektu)
- [Sposób działania](#sposób-działania)
- [Przykładowe użycie](#argumenty-do-wiersza-poleceń)
- [Plany rozwoju](#plany-rozwoju)
- [Współpraca i GitHub](#współpraca-i-github)

## Wymagania

- **Python 3.12.7** (lub wyższa)
- **MGLTools 1.5.7** – do przygotowania receptorów i ligandów
- **AutoDock4/AutoGrid4** oraz **AutoDock-GPU** (opcjonalnie)
- **MODELLER 10.6** 
- Biblioteki Python, takie jak:
  - `tqdm`
  - `openbabel-wheel 3.1.1.21`
  - `python-dotenv`
  - `pdbfixer`
  - `numpy`
  - `openmm`

## Instalacja

1. **Sklonuj repozytorium:**
   ```bash
   git clone https://github.com/siemania/mgr_AM_JL.git
   ```
2. **Przejdź do katalogu projektu:**
   ```bash
   cd mgr_AM_JL
   ```
3. **Utwórz środowisko wirtualne (opcjonalnie):**
   ```bash
   python -m venv .venv
   # Windows:
   .\.venv\Scripts\activate
   # Linux/Mac:
   source .venv/bin/activate
   ```
4. **Zainstaluj wymagane biblioteki:**
   ```bash
   pip install -r requirements.txt
   ```
5. **Zainstaluj Modeller:**
   1) Zarejestruj się na stronie, w celu uzyskania klucza licencyjnego
    [Modeller Registration](https://salilab.org/modeller/registration.html)
   2) Pobierz odpowiednią wersje Modellera
   [Download Modeller](https://salilab.org/modeller/download_installation.html)
   3) Uruchom instalator `.exe` (tu będzie potrzebny klucz licencyjny)
   4) Utwórz plik `license.py` oraz dodaj go do folderu modlib. Przykład co powinien zawierać taki plik:
    ```bash  
   license = 'XXXXXXXXXXXX'
    ```
   5) Skonfiguruj zmienne środowiskowe (jeżeli konieczne)
   
   Przykład:
   * Windows (cmd):
    ```bash
        set PYTHONPATH "C:\Program Files\Modeller10.6\modlib"
        set PATH "%PATH%;C:\Program Files\Modeller10.6\bin"
    ```
   * Linux:
    ```bash
        export PYTHONPATH="/opt/modeller10.6/modlib"
        export PATH="$PATH:/opt/modeller10.6/bin"
    ```
6. **Dla posiadaczy Linuxa:**
   1) **Zainstaluj MGLTools:**
   2) Przejdź do katalogu projektu, np. `docking`:
      ```bash
      cd ~/ścieżka/do/projektu/docking
      ```
   3) Pobierz instalator ze strony MGLTools:  
      [MGLTools Download (Linux)](https://ccsb.scripps.edu/mgltools/downloads/)  
      *(Wybierz odpowiednią wersję np. `mgltools_x86_64Linux2_1.5.7.tar.gz (Linux 64 tarball installer 108Mb)` lub inną)*
   4) Rozpakuj archiwum:
      ```bash
      tar -xvzf mgltools_*.tar.gz
      ```
   5) Przejdź do folderu rozpakowanego instalatora i uruchom instalację:
      ```bash
      cd mgltools_*/   # np. cd mgltools_x86_64Linux2_1.5.7/
      ./install.sh
      ```
   6) Po zakończeniu instalacji, zmień nazwę folderu:
      ```bash
      cd ..
      mv mgltools_*/ MGLTools_Linux/
      ```
      
## Uruchomienie programu

### Pobieranie i przygotowanie plików `.pdb`

1. **Przygotowanie listy PDB:**
   - W pliku `.csv` (np. `baza_ids.csv`) w katalogu głównym projektu, w kolumnie o nazwie `id` umieść kody PDB interesujących Cię kompleksów, każdy w osobnej linii.

2. **Pobieranie struktur:**
   - W katalogu głównym projektu uruchom skrypt `request_pdb_files.py`, który automatycznie pobierze wskazane struktury PDB.
   - Po zakończeniu działania skryptu upewnij się, że w folderze `Docking/pdb_files/` znajdują się pobrane pliki `.pdb` z bazy [RCSB](rcsb.org).

3. **Naprawa i standaryzacja plików:**
   - Aby przygotować pliki do dalszej analizy (np. usunąć ligandy, uzupełnić atomy wodoru, wyczyścić błędy strukturalne), w folderze `Docking/` uruchom skrypt `fixing_pdb_files.py`.
   - Skrypt utworzy nowy folder `Docking/fixed_pdb/` zawierający zminimalizowane i gotowe do dokowania struktury.

4. **Dokowanie:**
   - Uruchom skrypt `dock.py`.
   - Skrypt przyjmuje jako wejście:
     - pliki `.pdb`,
     - folder z wieloma plikami `.pdb` jak `pdb_files/`, lub
     - plik `.txt` zawierający czteroliterowe kody PDB (każdy w osobnej linii), np. `list_of_ids.txt`.

5. **(Opcjonalnie) Analiza danych BioLiP:**
   - Po pobraniu pliku `BioLiP.txt` (np. z oficjalnej [strony](https://zhanggroup.org/BioLiP/) bazy), możesz skorzystać ze skryptu
     `binding_energy_reader.py` w głównym katalogu do automatycznej analizy energii wiązania:
     
     ```bash
     python binding_energy_reader.py -c baza_ids.csv -b BioLiP.txt -o pdb_energy/ligands_energy.txt
     ```

6. **(Opcjonalnie) Zbieranie wartości Energii i RMSD:**
   - W folderze `statistics/` znajdują się użyteczne narzędzia do wykorzystania z użyciem CLI w odpowiedniej kolejności
   - Skrypt `statystyka.py` zawiera wszystkie kroki aby otrzymać wykresy RMSD oraz korelacyjny dla **2 zbiorów** wyników dokowania (plików z rozszerzeniem .dlg),
     należy tylko do 2 osobnych folderów pogrupować swoje pliki `.dlg`
   
   - `extract_experiment_energy.py` służy do zebrania tylko niezbędnych wartości energii kcal/mol po obróbce `binding_energy_reader.py` na podstawie własnych plików `.pdb`
   - `extract_energy_values.py` służy do zebrania wartości energii kcal/mol oraz RMSD z wyników dokowań `.dlg` w formie tabeli do pojedynczego pliku `.txt`
   - `rmsd_histogram.py` służy do wykonania histogramów RMSD na zbiorach dokowań i porównania wartości między sobą
   - `energy_correlation.py` służy do wykonania wykresów korelacyjnych dla plików po obróbkach: `extract_energy_values.py` dla zbiorów dokowań
     oraz `extract_experiment_energy.py` dla eksperymentalnych wartości energii kcal/mol

   Przykłady użycia:
   ```bash
   python extract_energy_values.py -d docking_results/ -o wyniki_dokowania.txt
   python extract_experiment_energy.py -f ligands_energy.txt -d results_pdb -o experiment_energy.txt
   python rmsd_hist.py -s wyniki_standard.txt -b wyniki_fixed.txt
   python correlation.py -s wyniki_standard.txt -b wyniki_fixed.txt -f wyniki_eksperymentalne_kcalmol.txt -o wynik_korelacji.png
   ```
   | Argument | | Description |
   | ------ | -- | ------------------- |
   | --help | -h | Wyświetlenie pomocy |
   | --file | -f | Pliki wynikowe .dlg lub .txt |
   | --names | -n | Nazwy zbiorów dla programu `statystyka.py` |
   | --directory | -d | Wybór katalogu z plikami .dlg lub .pdb |
   | --standard | -s | Plik .txt z wynikami po obróbce |
   | --fixed | -b | Plik .txt z wynikami po obróbce |
   | --labels | -l | Nazwy zbiorów dla dokowań |

> **Uwaga:** Upewnij się, że wszystkie ścieżki do katalogów i plików są poprawnie ustawione oraz że wszystkie wymagane zależności zostały wcześniej zainstalowane.

## Struktura projektu

```
mgr_AM_JL/
├── README.md                    # Ten plik
├── .gitignore                   # Plik ignorujący pliki wynikowe i lokalne konfiguracje
├── requirements.txt             # Lista zależności (m.in. tqdm, python-dotenv)
├── request_pdb_files.py         # Pobiera z bazy RCSB.org pliki pdb na podstawie osobnego pliku .csv z id kompleksów w kolumnie
├── binding_energy_reader.py     # Przelicza wartości eksperymentalne na kcal/mol z bazy BioLip.txt
└── Docking/                     # Folder główny projektu
    ├── dock.py                  # Główny skrypt wywołujący automatyzację dokowania
    ├── fixing_pdb_files.py      # Skrypt do ulepszenia plików PDB
    ├── pdb_files/               # Pliki PDB receptorów (wejściowe)
    ├── ligands/                 # Pliki ligandów (wejściowe)
    ├── pdbqt_files/             # Wynikowe pliki PDBQT (nie commitowane)
    ├── grid_dock_files/         # Pliki parametrów siatki i logi AutoGrid/AutoDock (nie commitowane)
    ├── output_files/            # Wynikowe kompleksy receptor-ligand (nie commitowane)
    ├── statistics/              # Folder z pomocnymi skryptami do statystyki automatycznej lub manualnej
    └── *.py                     # Mnóstwo innych skryptów (AutoDockTools oraz własnych)
```

## Sposób działania

Projekt automatycznie wykonuje następujące kroki:
- **Przygotowanie receptora:**  
  Używając MGLTools, receptor (plik PDB) jest przekształcany do formatu PDBQT, dodawane są ładunki oraz inne niezbędne modyfikacje.
- **Przygotowanie ligandu:**  
  Podobnie ligand jest przygotowywany do dokowania.
- **Generacja plików parametrów dokowania:**  
  Tworzony jest plik GPF, który jest zapisywany, a następnie wykorzystywany przez AutoGrid/AutoDock.
- **Dokowanie:**  
  Uruchamiane są odpowiednie polecenia AutoDock lub AutoDock-GPU.
- **Minimalizacja receptora:**
  Do poprawienia konformacji receptora (m.in. dodania wodoru i minimalizacji energii) wykorzystywane są wcześniej wspominane moduły. W ramach dalszego rozwoju chcemy dodać automatyczną optymalizację konformacyjną wybranych reszt (His, Gln, Asn) oraz innych aspektów.

Wszystkie etapy są wywoływane automatycznie, a postęp monitorowany jest przy użyciu `tqdm`. Po zakończeniu, wyniki (przygotowane pliki receptorów, ligandów oraz kompleksów) są zapisywane w odpowiednich folderach, które nie są wersjonowane.

## Sposób działania - fixing_pdb_files.py

Ten skrypt automatyzuje proces uzupełniania i optymalizacji struktur białek zapisanych w plikach PDB. Korzysta z biblioteki Modeller oraz narzędzia Open Babel do dodawania atomów wodoru i przeprowadzania optymalizacji geometrii cząsteczek.

**Główne kroki działania:**
1. Przygotowanie środowiska i struktur
2. Uzupełnienia brakujących fragmentów w strukturze
3. Dodanie atomów wodoru
4. Poprawa orientacji niektórych reszt aminokwasowych
5. Optymalizacja geometrii

**Kluczowe klasy i funkcje**
+ `MyModel` - klasa rozszerzająca `AutoModel` z Modellera, pozwalająca na póżniejszy wybór atomów do modelowania
+ `PDBModelOptimization` - klasa zarządzająca całym procesem
  + `prepare_alignment` - tworzy dopasowanie sekwencji i struktury
  + `cleanup_working_files` - usuwa na bieżąco pliki folderu roboczego oprócz alignment.
  + `add_hygrodens` - dodaje atomy wodoru za pomocą OpenBabel
  + `flip_residues` - koryguje orientacje reszt ASN, GLN i HIS
  + `optimize_heavy_atom` - optymalizuje geometrię atomów ciężkich
  + `optimize_full_structure` - pełna optymalizacja struktury (gradienty + dynamika molekularna)
  + `fill_missing_residues_and_atoms` - główna metoda przetwarzająca wszystkie pliki PDB z folderu wejściowego

## Argumenty do wiersza poleceń

  ```bash
  python dock.py [opt. arguments]
  python dock.py -f 1akw.pdb 1hgw.pdb
  python dock.py -d pdb_files/ -s "receptor" "ligand" "grid" "autogrid" "dock"
  ```
   
| Argument | | Description | Default value |
| ------ | -- | ------------------- | -- |
| --help | -h | Wyświetlenie pomocy | No |
| --file | -f | Pliki kompleksów w formacie .pdb | each pdb_file in pdb_files/ |
| --directory | -d | Wybór katalogu z plikami .pdb | pdb_files/ |
| --select_command | -s | Wybór komend do wykonania:<br> receptor, ligand, grid, autogrid, dock, autodock, complex | All |

## Plany rozwoju

W kolejnych etapach projektu planujemy rozwój modułów optymalizacyjnych, takich jak:
- Opcjonalne dokowanie z wodą lub bez wody.
- Dodanie interfejsu graficznego lub narzędzi do analizy wyników.

## Współpraca i GitHub
Program został stworzony na potrzeby pracy magisterskiej pod tytułem:  
> Procedura automatycznego przygotowania receptorów do symulacji dokowania  

Uczelnia macierzysta: [Politechnika Gdańska](https://pg.edu.pl)

Za program odpowiedzialni są:  
[inż. Jan Latt](https://www.github.com/Tajgero)  
[inż. Anna Muńko](https://github.com/siemania)

Projekt jest hostowany na GitHubie i współpracujemy przy użyciu PyCharm:
- **.gitignore:** Wszystkie pliki wynikowe, lokalne konfiguracje oraz foldery generowane automatycznie (takie jak `pdbqt_files/`, `output_files/`, `grid_dock_files/`, `ligands/`, `work_folder/`, `fixed_pdb`) są ignorowane przez Git.
