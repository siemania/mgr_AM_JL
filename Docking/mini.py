import os
from modeller import *
from modeller.optimizers import ConjugateGradients, MolecularDynamics, actions
from modeller.scripts import complete_pdb



def fill_missing_atoms(input_dir, output_dir):

    log.level(output=1, notes=1, warnings=1, errors=1, memory=0)

    # Ustawianie srodowiska
    env = Environ()

    env.io.atom_files_directory = input_dir
    env.edat.dynamic_sphere = True

    env.libs.topology.read('${LIB}/top_allh.lib')
    env.libs.parameters.read('${LIB}/par.lib')


    for filename in os.listdir(input_dir):
        if not filename.endswith(".pdb"):
            continue

        name = os.path.splitext(filename)[0]
        code = filename
        template = os.path.join(PROJECT_ROOT, "ent_files", name + ".entOrg")

        if not os.path.isfile(template):
            print(f"Brak pliku szablonu (entOrg) dla {name}, pomijam...")
            continue

        print(f"Przetwarzanie {filename}...")

        try:

            #Wczytanie modelu bedacego szablonem
            model_T = complete_pdb(env, template)
            chain_T = model_T.chains[0]
            resids_T = set(res.num for res in chain_T.residues)

            #Uzupelnienie modelu do poprawy i nalozenie restrainow
            model = complete_pdb(env, code)
            chain = model.chains[0]
            resids = set(res.num for res in chain.residues)

            missing_resids = resids_T - resids

            if not missing_resids:
                print(f"Brak brakujących reszt w {name}, pomijam.")
                continue

            select_T = Selection(*[res for res in chain_T.residues if res.num in missing_resids])
            select = Selection(*[res for res in chain.residues if res.num in missing_resids])

            #Restreiny szablonu dla brakujacych reszt
            restrains_T = model_T.restraints
            restrains_T.make(select_T, restraint_type='dihedral', spline_on_site=False)
            restrains_T.write(file=code + '.rsrT')


            #Restreiny modelu dla brakujacych reszt
            restrains = model.restraints
            restrains.make(Selection(model), restraint_type='stereo', spline_on_site=False)
            restrains.write(file=code + '.rsr1')
            restrains.append(file=code + '.rsrT')
            restrains.write(file=code + '.rsr2')

            #Optymalizacja
            cg = ConjugateGradients(output='REPORT')
            md = MolecularDynamics(output='REPORT')

            #Utworzenie pliku prostego loga
            log_path = os.path.join(output_dir, name + '.D00000001')
            with open(log_path, 'w') as trcfil:

                cg.optimize(select, max_iterations=20, actions=actions.Trace(5, trcfil))

                md.optimize(select,
                        temperature=300, max_iterations=50,
                        actions=[actions.WriteStructure(10, code + '.D9999%04d.pdb'),
                                actions.Trace(10, trcfil)])

                cg.optimize(select, max_iterations=20, actions=[actions.Trace(5, trcfil)])

            #Zapis koncowego modelu
            model.write(file=os.path.join(output_dir, code + 'Bcg2'))
            print(f"Zapisano do: {output_dir}")

        except Exception as e:
            print(f"Błąd podczas przetwarzania {filename}: {e}")


if __name__ == "__main__":
    username = os.environ.get('USERNAME') or os.environ.get('USER')
    PROJECT_ROOT = rf"C:\Users\{username}\PycharmProjects\mgr_AM_JL"
    INPUT_DIR = os.path.join(PROJECT_ROOT, 'pdb_files')
    OUTPUT_DIR = os.path.join(PROJECT_ROOT, 'fixed_pdb')

    os.makedirs(OUTPUT_DIR, exist_ok=True)


    fill_missing_atoms(INPUT_DIR, OUTPUT_DIR)

