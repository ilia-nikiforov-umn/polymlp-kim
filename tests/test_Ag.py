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
        "Ag4",
        cell=axis,
        scaled_positions=scaled_positions,
        pbc=True,
    )
    return atoms


def get_structure2():
    """Return test structure."""
    axis = np.array(
        [
            [4.05, 0.01, -0.01],
            [4.05, 4.06, 0.02],
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


def test_eval2():
    """Test property calculations."""
    model = "Polymlp_Seko_2024_Ag__MO_010001100136_000"
    atoms = get_structure2()
    energy, forces, stress = eval_using_kim(model, atoms)
    e_true = 4.198900014523197
    f_true = [[ 2.04197992, -8.38688164,  0.04647052],
              [-0.83417846,  4.26332429,  0.10116097],
              [ 0.97160298, -4.3146287 , -0.07354999],
              [-2.17940444,  8.43818605, -0.07408151]]
    s_true = [
        5.49418388e+00,  
        3.06901857e+02,  
        5.21703045e+00, 
        -2.90077394e+00,
        2.22746389e+00,  
        3.06713742e-02,
    ]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval3():
    """Test property calculations."""
    model = "Polymlp_Seko_2024p1hybrid_Ag__MO_000000000000_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)
    print(energy)
    print(forces)
    print(stress)
    e_true = -9.990874169902824
    f_true = [[-0.241219  ,-0.52440161,-0.05992873],
              [-0.00514723, 0.40081802, 0.02103461],
              [ 0.11136887,-0.28496934, 0.0589867 ],
              [ 0.13499736, 0.40855293,-0.02009258]]
    s_true = [
        7.39887587e+00, 
        7.29192330e+00, 
        7.20930295e+00, 
        -1.59844219e-03,
         -2.57000240e-01,
         1.56510639e-01,
    ]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)
