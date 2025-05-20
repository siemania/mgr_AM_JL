# Automatyzacja Dokowania

Projekt **mgr_AM_JL** to zestaw skryptów mających na celu automatyzację przygotowania receptorów i ligandów do dokowania. Początkowo projekt automatycznie wykonuje serię poleceń związanych z przygotowaniem receptorów i ligandów przy użyciu MGLTools oraz AutoDock/AutoDock GPU, a w dalszej perspektywie planowane jest rozwinięcie narzędzi do:

- Optymalizacji konformacji wybranych reszt (np. HIS, GLU, ASP)
- Dokowania z wodą i bez wody
- Optymalizacji układu wodorów (obrót grup metylowych itp.)
- Dodatkowych algorytmów optymalizacyjnych wykorzystujących pola siłowe z programów jak np. OpenMM, OpenBabel, GROMACS

## Spis treści

- [Wymagania](#wymagania)
- [Instalacja](#instalacja)
- [Konfiguracja](#konfiguracja)
- [Struktura projektu](#struktura-projektu)
- [Sposób działania](#sposób-działania)
- [Przykładowe użycie](#przykładowe-użycie)
- [Plany rozwoju](#plany-rozwoju)
- [Współpraca i GitHub](#współpraca-i-github)

## Wymagania

- **Python 3.12.7** (lub wyższa)
- **MGLTools 1.5.7** – do przygotowania receptorów i ligandów
- **AutoDock4/AutoGrid4** oraz **AutoDock-GPU** (opcjonalnie)
- Biblioteki Python, takie jak:
  - `tqdm`
  - `openbabel-wheel`
  - `python-dotenv`
  - Standardowe biblioteki: `os`, `glob`, `time`, `datetime`, `subprocess`, `sys`

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
<!--   ALBO ZMODYFIKOWAĆ ALBO ZMIENIĆ, BO NIE WIADOMO CZY BĘDZIE POTRZEBNE
## Konfiguracja

- **Plik `.env`:**  
  Utwórz lokalny plik `.env` (który nie jest commitowany) na podstawie szablonu `.env.example`. Plik ten powinien zawierać zmienne środowiskowe specyficzne dla Twojego systemu, np. ścieżki do interpretera Pythona, AutoGrid, AutoDock czy AutoDock-GPU.  
  Przykład:
  ```dotenv
  PYTHON_PATH=MGLTools/python.exe
  AUTOGRID_PATH=autogrid4.exe
  AUTODOCK_PATH=autodock4.exe
  AUTODOCK_GPU_PATH=./bin/autodock_gpu_128wi
  ```
-->
## Struktura projektu

```
mgr_AM_JL/
├── README.md              # Ten plik
├── .gitignore             # Plik ignorujący pliki wynikowe i lokalne konfiguracje
├── requirements.txt       # Lista zależności (m.in. tqdm, python-dotenv)
├── .env.example           # Szablon zmiennych środowiskowych
└── Docking/               # Folder główny projektu
    ├── _INFO.txt              # Informacje o dock.py (część dokumentacyjna prywatne zmiany)
    ├── TEST_SKRYPTU/          # Folder służący do testowania funkcji na kilku plikach
    ├── dock.py                # Główny skrypt wywołujący automatyzację dokowania
    ├── pdb_files/             # Pliki PDB receptorów (wejściowe)
    ├── ligands/               # Pliki ligandów (wejściowe)
    ├── pdbqt_files/           # Wynikowe pliki PDBQT (nie commitowane)
    ├── grid_dock_files/       # Pliki parametrów siatki i logi AutoGrid/AutoDock (nie commitowane)
    ├── output_files/          # Wynikowe kompleksy receptor-ligand (nie commitowane)
    └── modify_parameters.py   # Moduł modyfikujący parametry (np. ścieżki, automatyzację)
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
  Do poprawienia konformacji receptora (m.in. dodania wodoru i minimalizacji energii) wykorzystywane są wcześniej wspominane moduły. W ramach dalszego rozwoju chcemy dodać automatyczną optymalizację konformacyjną wybranych reszt (HIS, GLU, ASP) oraz innych aspektów.

Wszystkie etapy są wywoływane automatycznie, a postęp monitorowany jest przy użyciu `tqdm`. Po zakończeniu, wyniki (przygotowane pliki receptorów, ligandów oraz kompleksów) są zapisywane w odpowiednich folderach, które nie są wersjonowane.

## Przykładowe użycie

1. **Skonfiguruj środowisko** (utwórz `.env`, zainstaluj zależności).
2. **Uruchom główny skrypt:**
   ```bash
   python dock.py
   ```
3. **Monitoruj postęp:**  
   Pasek postępu pojawi się w terminalu, a komunikaty zostaną wyświetlone na bieżąco.

## Dodatkowe argumenty
| Argument | | Description | Default value |
| ----------- | - | ----------- | ---------- |
| --file | -f | Pliki kompleksów w formacie .pdb |

## Plany rozwoju

W kolejnych etapach projektu planujemy:
- Rozwój modułów optymalizacyjnych, takich jak:
  - Minimalizacja konformacji receptorów za pomocą OpenBabel, RDKit, GROMACS.
  - Automatyczne ustawianie wiązań (constraints) na ciężkie atomy, aby ich pozycje były modyfikowane łagodniej.
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
- **.gitignore:** Wszystkie pliki wynikowe, lokalne konfiguracje oraz foldery generowane automatycznie (takie jak `pdbqt_files/`, `output_files/`, `grid_dock_files/`, `ligands/`) są ignorowane przez Git.
