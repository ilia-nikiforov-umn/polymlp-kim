"""Tests of property calculations in Ag-Au."""

import numpy as np

from ase.atoms import Atoms
from eval_kim import eval_using_kim


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
        "CuAgAu2",
        cell=axis,
        scaled_positions=scaled_positions,
        pbc=True,
    )
    return atoms


def test_eval1():
    """Test property calculations."""
    model = "Polymlp_Seko_2024_CuAgAu__MO_030000020182_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -12.179284644330645
    f_true = [[-0.09204175, -0.23351199, -0.02565584],
              [-0.00185949,  0.28817092,  0.01279268],
              [ 0.0393424 , -0.35178773,  0.03880553],
              [ 0.05455884,  0.2971288 , -0.02594236]]
    s_true = [
        -1.3434409, -1.31555145, -1.40991519, -0.00826187, -0.18569514, 0.10744792
    ]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)
