# Uzupelnia plik PDB z dziurami o brakujace aminokwasy
# w tej wersji nie rusza pozostalych aminokwasow, ktorych
# pozycje sa znane

import sys
from modeller import *
from modeller.scripts import complete_pdb

filein  = sys.argv[1]
fileout = filein + '.fixed.pdb'

log.verbose()
env = environ()

# Uzyj topologii zawierajacej wodory
env.libs.topology.read('${LIB}/top_allh.lib')
env.libs.parameters.read('${LIB}/par.lib')

env.io.water = True

# Uzupelnij brakujace atomy, w tym wodory
mdl = complete_pdb(env, filein)

# Zapisz wynik
mdl.write(file=fileout)
