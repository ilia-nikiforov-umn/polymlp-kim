"""Compute the cohesive energy, forces, and pressure of given POSCAR."""

import numpy as np

from ase.calculators.kim import KIM
from ase.atoms import Atoms
from ase.units import GPa



models = [
    "Polymlp_Seko_2022_AgAu__MO_020000020182_000",
	"Polymlp_Seko_2022_AgAu__MO_020000120203_000",
	"Polymlp_Seko_2022_AgAu__MO_020000220257_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000320269_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000420416_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000520441_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000620722_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000720725_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000820734_000",
#	"Polymlp_Seko_2022_AgAu__MO_020000910050_000",
#	"Polymlp_Seko_2022_AgAu__MO_020001010051_000",
]

axis = np.array(
    [
        [4.05, 0.00, 0.00],
        [0.00, 4.05, 0.00],
        [0.00, 0.00, 4.05],
    ]
)
scaled_positions = np.array(
    [
        [0.01,0.02,0.002],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]
)
atoms = Atoms(
    "Ag2Au2",
    cell=axis,
    scaled_positions=scaled_positions,
    pbc=True,
)

for model in models:
    calc = KIM(model)
    atoms.calc = calc

    ecoh = -atoms.get_potential_energy() 
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    pressure_GPa = (-sum(stress[:3]) / 3.0) / GPa

    print("----------------", model, "----------------")
    print("Computed cohesive energy of {:.12f} eV/atom".format(ecoh))
    print("Computed pressure of {:.9f} GPa".format(pressure_GPa))
    print("Computed forces")
    print(forces)
    print("Computed stress")
    print(-stress / GPa)
    print("")
