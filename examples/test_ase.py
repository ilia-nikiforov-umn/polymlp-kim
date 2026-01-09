#!/usr/bin/env python3
"""
Compute the cohesive energy and pressure of an FCC Al crystal using the
Ercolessi-Adams EAM potential implemented as a Portable Model (PM) in
OpenKIM for the experimental lattice constant a0=4.05 Angstrom.
"""
from ase.calculators.kim import KIM
from ase.lattice.cubic import FaceCenteredCubic
from ase.units import GPa

# Set up crystal and calculator
a0 = 4.05  # experimental lattice constant
atoms = FaceCenteredCubic("Ag", latticeconstant=a0)
#calc = KIM("Tersoff_LAMMPS_Tersoff_1989_SiC__MO_171585019474_004")
calc = KIM("Polymlp_Seko_2022_AgAu__MO_000000111111_000")
atoms.set_calculator(calc)

# Compute energy/pressure
# ecoh = -atoms.get_potential_energy() / len(atoms)
ecoh = -atoms.get_potential_energy() 
stress = atoms.get_stress()
pressure_MPa = (-sum(stress[:3]) / 3.0) * 1e3 / GPa

print("Computed cohesive energy of {:.10f} eV/atom".format(ecoh))
print("Computed pressure of {} MPa".format(pressure_MPa))
