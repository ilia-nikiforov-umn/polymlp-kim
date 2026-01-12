"""Functions for evalating properties using KIM and ASE."""

import numpy as np

from ase.calculators.kim import KIM
from ase.atoms import Atoms
from ase.units import GPa


def eval_using_kim(model: str, atoms: Atoms):
    """Evaluate properties using KIM."""
    calc = KIM(model)
    atoms.calc = calc

    ecoh = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    pressure_GPa = (-sum(stress[:3]) / 3.0) / GPa

    # Order (xx, yy, zz, yz, xz, xy) in ASE -> (xx, yy, zz, xy, yz, zx) in pypolymlp.
    pypolymlp_stress = - stress[np.array([0, 1, 2, 5, 3, 4])] / GPa
    return ecoh, forces, pypolymlp_stress
