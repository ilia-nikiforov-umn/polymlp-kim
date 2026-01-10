#!/usr/bin/env python3
"""
Compute the cohesive energy and pressure of an FCC Al crystal using the
Ercolessi-Adams EAM potential implemented as a Portable Model (PM) in
OpenKIM for the experimental lattice constant a0=4.05 Angstrom.
"""
from ase.calculators.kim import KIM
from ase.lattice.cubic import FaceCenteredCubic
from ase.units import GPa

models = [
    "Polymlp_Seko_2022_AgAu__MO_000000000001_000",
    "Polymlp_Seko_2022_AgAu__MO_000000111111_000",
    "Polymlp_Seko_2022_AgAu__MO_000000000325_000",
]

a0 = 4.05  # experimental lattice constant
atoms = FaceCenteredCubic("Ag", latticeconstant=a0)
for model in models:
    calc = KIM(model)
    atoms.calc = calc

    ecoh = -atoms.get_potential_energy() 
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    pressure_GPa = (-sum(stress[:3]) / 3.0) / GPa

    print("######", model, "#####")
    print("Computed cohesive energy of {:.12f} eV/atom".format(ecoh))
    print("Computed pressure of {:.9f} GPa".format(pressure_GPa))
    print("Computed forces")
    print(forces)
    print("Computed stress")
    print(-stress / GPa)
