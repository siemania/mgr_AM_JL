# przykladowa optymalizacja fragmentu struktury z dodatkowymi wiezami
# wzietymi z innej struktury

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import ConjugateGradients, MolecularDynamics, actions


log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
env = Environ()
env.io.atom_files_directory = ['.']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')


####################################
# nazwy plikow
####################################
code        = 'pdb1crn' # ten bedzie poprawiany. Recznie wykasowalem w nim amino 33 i 44 
template    = 'pdb1crn.entOrg'

####################################
###  Zapis restreinow z na podstawie struktury "szablonu" 
####################################
mdl_T  = complete_pdb(env, template)
rsr_T  = mdl_T.restraints

# tylko aminokwasy 33 i 44 beda optymalizowane. Numeracja jako str (czyli '33') oznacz
# numeracje jak w pliku pdb. Numeracja jako int (czyli 33) oznacza numeracje pythonowa
# czyli pierwszy aminokwas ma numer 0.
sel_T  = Selection(mdl_T.chains[0].residues['44'], mdl_T.chains[0].residues['33'])
rsr_T.make(sel_T, restraint_type='dihedral', spline_on_site=False)
rsr_T.write(file=code+'.rsrT')

###
### Wczytanie modelowanego pliku
###
mdl = complete_pdb(env, code) # rozszerzenie dodawane jest automatycznie
aln = alignment(env)
rsr = mdl.restraints

atmsel  = Selection(mdl) # to w sumie nie jest potrzebne. Jest tylko jako przyklad
ressel  = mdl.chains[0].residues
sele    = Selection(ressel['44'], ressel['33'])


rsr.make(atmsel, restraint_type='stereo', spline_on_site=False)
# powyzej sa restreiny dla wszystkich atomow (atmsel) ale w sumie 
# wystarczylyby tylko dla zaznaczonych 33 i 44 (sele)
# zapisane sa dla porownania restreiny przerd i po dodaniu tych z szablonu
rsr.write(file=code+'.rsr1')
rsr.append(file=code+'.rsrT')
rsr.write(file=code+'.rsr2')

# teraz bedzie optymalizacha wybranych aminokwasow na kilka sposobow jako przyklady
cg = ConjugateGradients(output='REPORT')
md = MolecularDynamics(output='REPORT')

# Utworzenie pliku prostego loga
trcfil = open(code+'.D00000001', 'w')

# Minimalizacja CG dla zaznaczonych atomow; Zapis co 5 krokow
cg.optimize(sele, max_iterations=20, actions=actions.Trace(5, trcfil))
mdl.write(file=code+'.Bcg')

# Dla przykladu dynamika z zapisem kilku klatek co 10 krokow ('1fas.D9999xxxx.pdb')
md.optimize(sele, 
            temperature=300, max_iterations=50,
            actions=[actions.WriteStructure(10, code+'.D9999%04d.pdb'),
                     actions.Trace(10, trcfil)])

# Jeszcze jedna optymalizacja ostatniej klatki po dynamice
cg.optimize(sele, 
            max_iterations=20,
            actions=[actions.Trace(5, trcfil)])

mpdf = atmsel.energy()

mdl.write(file=code+'.Bcg2')

