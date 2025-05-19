import numpy as np

def calculate_offset(reference, target, gpf_file):
    """
    Oblicza przesunięcie (offset) między koordynatami docelowymi i referencyjnymi
    używając NumPy dla wydajnych obliczeń wektorowych.
    """
    # Settings
    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
    
    # Wyodrębnij tylko koordynaty XYZ z referencji (ignoruj ID w ostatniej kolumnie)
    ref_coords = np.array([r[:3] for r in reference])
    target_coords = np.array(target)
    gpf_coords = np.array(gpf_file)
    
    # Oblicz offset jako target - reference
    offsets = target_coords - ref_coords
    
    # Dodaj offset do wiersza "gridcenter" grida
    gridcenter = gpf_coords + offsets
    
    # Wyświetl wyniki dla każdej pary
    print("Offsety dla wszystkich struktur:")
    for number, (ref, tgt, offset, gpf, grid) in enumerate(zip(reference, target, offsets, gpf_file, gridcenter)):
        print(f"Para {number+1} ({ref[3]}):")
        print(f"  Reference: {ref[:3]}")
        print(f"  Target:    {tuple(tgt)}")
        print(f"  Offset:    {offset}")
        print(f"  Gpfcoords: {gpf}")
        print(f"  Grid:      {grid}\n")

    # # Odczytaj zawartość pliku
    # with open(gpf_file, 'r') as file:
    #     lines = file.readlines()

if __name__ == "__main__":
    
    reference = [(34.149, -13.537, -32.812, "1trd"),
                 (-12.427,  87.778,  34.508, "1vot"),
                 (18.199,  21.058,  70.786, "2fzz"),
                 (3.353,  -4.178,  -2.630, "2jfh"),
                 (11.040,  22.042,  10.669, "4eu3"),
                 (42.255, -11.646,  34.367, "4wxi"),
                 (7.041, -13.169, -29.389, "5diu"),
                 (-5.519,  -0.931,  20.437, "5i3h"),
                 (-22.086,  -1.720, -24.347, "5j7q"),
                 (19.561, -16.549,  -2.236, "5uex")]
    
    target = [(67.960,  43.660,  53.820),
              (24.880,  67.110,  19.840),
              (53.780,  34.080,  36.380),
              (19.300,  46.780,  22.370),
              (78.620,  75.630, 100.500),
              (43.050,  40.930,  33.670),
              (33.490,  15.150,  15.030),
              (27.440,  35.180,  84.520),
              (38.900,  29.060,  22.690),
              (28.280,  13.880,  15.930)]
    
    gpf = [(7.195, -23.470, -17.116),
           (3.186, 68.014, 63.086),
           (21.349, 30.309, 61.055),
           (20.993, -1.904, 20.391),
           (-20.260, -25.308, -29.839),
           (41.728, -20.058, 29.363),
           (12.923, -3.459, -7.720),
           (25.562, -5.893, -24.569),
           (-18.697, -0.573, -24.975),
           (12.218, 12.346, 19.473)]
    
    calculate_offset(reference, target, gpf)
    
    
    
    