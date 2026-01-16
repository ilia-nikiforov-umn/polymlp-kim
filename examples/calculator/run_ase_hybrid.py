"""Compute the cohesive energy, forces, and pressure of given POSCAR."""

import numpy as np

from ase.calculators.kim import KIM
from ase.atoms import Atoms
from ase.units import GPa
from ase.build.supercells import make_supercell


models = [
    "Polymlp_Seko_2024p1hybrid_Ag__MO_000000000000_000"
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
        [0.010, 0.020, 0.002],
        [0.000, 0.500, 0.500],
        [0.501, 0.010, 0.499],
        [0.500, 0.500, 0.000],
    ]
)
atoms = Atoms(
    "Ag4",
    cell=axis,
    scaled_positions=scaled_positions,
    pbc=True,
)
atoms = make_supercell(atoms, np.eye(3, dtype=int) * 10)

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
