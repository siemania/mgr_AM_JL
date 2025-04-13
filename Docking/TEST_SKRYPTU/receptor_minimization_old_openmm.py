#!/usr/bin/env python
# -*- coding: utf-8 -*-

from openmm.app import *
from openmm import *
from openmm.unit import *
import os
import sys


class ReceptorMinimizerOpenMM:
    def __init__(self, pdb_file):
        """
        Inicjalizacja – wczytanie receptora z pliku PDB (bez atomów wodoru).
        """
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"Plik {pdb_file} nie został znaleziony!")

        # Wczytujemy receptor; PDBFile domyślnie nie dodaje wodoru
        self.pdb = PDBFile(pdb_file)
        self.modeller = Modeller(self.pdb.topology, self.pdb.getPositions())

    def add_hydrogens(self, pH=7.0):
        """
        Dodaje brakujące atomy wodoru do receptora, biorąc pod uwagę pH.
        """
        # Dodanie wodoru w oparciu o podane pH
        self.modeller.addHydrogens(pH=pH)
        print("Dodano atomy wodoru.", file=sys.stdout, flush=True)

    def minimize_with_constraints(self, force_constant=1000.0, tol=10 * kilojoule_per_mole):
        """
        Minimalizuje receptor, nakładając restrykcje na ciężkie atomy.

        Parametry:
         - force_constant: stała harmoniczna restrykcji (w kilodżulach/(mol*nm^2))
         - tol: tolerancja minimalizacji (np. 10 kJ/mol)
         - max_iters: maksymalna liczba iteracji (używana domyślnie przez LocalEnergyMinimizer)
        """
        # Używamy force fielda Amber14
        forcefield = ForceField("amber14-all.xml")
        # Tworzymy system bez domyślnych restrykcji
        system = forcefield.createSystem(self.modeller.topology,
                                         nonbondedMethod=NoCutoff,
                                         constraints=None)

        # Tworzymy custom force dla restrykcji: harmoniczny potencjał, który karze odchylenie od pozycji początkowych
        restraint_force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint_force.addGlobalParameter("k", force_constant * kilojoule_per_mole / nanometer ** 2)
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")

        positions = self.modeller.getPositions()

        # Dodajemy restrykcje tylko dla ciężkich atomów na końcach N i C (bez atomów wodoru)
        for chain in self.modeller.topology.chains():
            n_terminal_res = next(chain.residues())  # Pierwszy reszt - N-końcowy
            c_terminal_res = list(chain.residues())[-1]  # Ostatni reszt - C-końcowy

            for atom in n_terminal_res.atoms():
                if atom.element is not None and atom.element.symbol != "H":
                    idx = atom.index
                    pos = positions[idx]
                    restraint_force.addParticle(idx, [pos.x, pos.y, pos.z])

            for atom in c_terminal_res.atoms():
                if atom.element is not None and atom.element.symbol != "H":
                    idx = atom.index
                    pos = positions[idx]
                    restraint_force.addParticle(idx, [pos.x, pos.y, pos.z])

        system.addForce(restraint_force)

        # Konfiguracja minimalizacji przy użyciu integratora Verlet (używamy niewielkiego kroku, bo minimalizacja statyczna)
        integrator = VerletIntegrator(1.0 * femtosecond)
        platform = Platform.getPlatformByName("Reference")
        simulation = Simulation(self.modeller.topology, system, integrator, platform)
        simulation.context.setPositions(positions)

        print("Rozpoczynam minimalizację receptora z restrykcjami...", file=sys.stdout, flush=True)
        # Używamy LocalEnergyMinimizer; max_iters nie jest bezpośrednio przekazywany, ale możemy zignorować,
        # bo OpenMM używa wewnętrznych kryteriów
        LocalEnergyMinimizer.minimize(simulation.context, tolerance=tol)
        print("Minimalizacja zakończona.", file=sys.stdout, flush=True)

        # Aktualizujemy pozycje w modellerze
        state = simulation.context.getState(getPositions=True)
        self.modeller.positions = state.getPositions()

    def write_pdb(self, output_file):
        """
        Zapisuje zminimalizowany receptor do pliku PDB.
        """
        with open(output_file, "w") as f:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, f)
        print(f"Zminimalizowany receptor zapisany jako {output_file}", file=sys.stdout, flush=True)


if __name__ == '__main__':
    input_pdb = "pdbqt_files\\1hsg_receptor.pdb"
    output_pdb = "output_files\\1hsg_receptor_minimized.pdb"

    os.makedirs("output_files", exist_ok=True)

    # Utwórz obiekt minimalizatora i wykonaj kolejne kroki
    minimizer = ReceptorMinimizerOpenMM(input_pdb)
    minimizer.add_hydrogens(pH=7.0)
    minimizer.minimize_with_constraints(force_constant=1000.0, tol=10 * kilojoule_per_mole)
    minimizer.write_pdb(output_pdb)
