"""Tests of property calculations in Ag-Au."""

from pathlib import Path

import numpy as np

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
    return atoms


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


def test_eval1():
    """Test property calculations."""
    model = "Polymlp_Seko_2024_Ag__MO_010001100136_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -9.99070955606906
    f_true = [[-0.24055861, -0.52328104, -0.05976086],
              [-0.00554532,  0.40021048,  0.02094011],
              [ 0.11133271, -0.28492353,  0.05898258],
              [ 0.13477122,  0.4079941 , -0.02016183]]
    s_true = [
        7.43149206,  
        7.32045211, 
        7.23729695, 
        -4.91592175e-04, 
        -2.56989263e-01,  
        1.56227642e-01,
    ]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)
