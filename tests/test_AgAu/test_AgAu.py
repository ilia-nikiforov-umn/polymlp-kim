"""Tests of neighbor calculations."""

from pathlib import Path

import numpy as np
import pytest

from ase.calculators.kim import KIM
from ase.atoms import Atoms
from ase.units import GPa


cwd = Path(__file__).parent


def get_structure1():
    """Return test structure."""
    axis = np.array(
        [
            [4.05, 0.01, -0.01],
            [0.00, 4.06, 0.02],
            [0.00, 0.00, 4.07],
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
    return atoms


def test_eval1():
    """Test property calculations."""
    calc = KIM(model)
    atoms.calc = calc

    ecoh = -atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    pressure_GPa = (-sum(stress[:3]) / 3.0) / GPa

    pypolymlp_stress = - stress[np.array([0, 1, 2, 5, 4, 3])]
#     model = "Polymlp_Seko_2022_AgAu__MO_020000020182_000"
#     st = Poscar(cwd / "POSCAR").structure
#     prop = Properties(pot=cwd / model / "polymlp.yaml")
#     energy, forces, stresses = prop.eval(st)

    print(energy)
    print(forces)
    print(stresses)
    assert 1 == 0
