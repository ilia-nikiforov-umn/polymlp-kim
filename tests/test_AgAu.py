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
        "Ag2Au2",
        cell=axis,
        scaled_positions=scaled_positions,
        pbc=True,
    )
    return atoms


def test_eval1():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000020182_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.234256372894958
    f_true = [[-0.27460731, -0.57769267, -0.06669096],
              [-0.01834291,  0.44498066,  0.02475438],
              [ 0.13302284, -0.38203952,  0.07262356],
              [ 0.15992738,  0.51475153, -0.03068698]]
    s_true = [9.55151215, 8.75862221, 8.64869267, 0.03470536, -0.23675401, 0.15981854]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval2():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000120203_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.231262011084647
    f_true = [[-0.2755677 , -0.57783952, -0.06698578],
              [-0.01864271,  0.4441513 ,  0.02466768],
              [ 0.1332546 , -0.37779368,  0.07230968],
              [ 0.16095581,  0.5114819 , -0.02999158]]
    s_true = [10.02353709, 8.854922, 8.74257165, 0.0325456, -0.22551926, 0.15727892]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval3():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000220257_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.232157286398184
    f_true = [[-0.27298747, -0.58018695, -0.06715227],
              [-0.02033261,  0.44611096,  0.02481966],
              [ 0.13317646, -0.38362256,  0.07302424],
              [ 0.16014362,  0.51769855, -0.03069162]]
    s_true = [9.71383737, 8.93589813, 8.82995387, 0.03558284, -0.23397547, 0.15604278]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval4():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000320269_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.23258732954334
    f_true = [[-0.27338815, -0.57902214, -0.06685843],
              [-0.02068024,  0.44664909,  0.02502294],
              [ 0.13351985, -0.38351019,  0.07261065],
              [ 0.16054855,  0.51588324, -0.03077515]]
    s_true = [9.74056784, 8.96833653, 8.86361243, 0.03373103, -0.22822689, 0.15660628]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval5():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000420416_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.23548602102875
    f_true = [[-0.27353559, -0.57999621, -0.06714526],
              [-0.02002739,  0.44598892,  0.02485334],
              [ 0.13316029, -0.38210263,  0.07289718],
              [ 0.16040269,  0.51610992, -0.03060526]]
    s_true = [9.71460123, 8.94009381, 8.82934187, 0.03567645, -0.22755002, 0.15637431]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval6():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000520441_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.236685599020463
    f_true = [[-0.27448732, -0.57948433, -0.06706481],
              [-0.01957311,  0.44638108,  0.02503769],
              [ 0.13333192, -0.38075219,  0.07255791],
              [ 0.16072851,  0.51385544, -0.03053079]]
    s_true = [9.85958667, 8.99648287, 8.88426581, 0.02721249, -0.2414297, 0.16289111]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval7():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000620722_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.235391426480263
    f_true = [[-0.28002057, -0.58326896, -0.06785937],
              [-0.01729502,  0.44649802,  0.02474336],
              [ 0.13503931, -0.38140713,  0.07326037],
              [ 0.16227627,  0.51817807, -0.03014436]]
    s_true = [9.04246599, 8.60650644, 8.49134125, 0.03140227, -0.23394861, 0.16046488]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval8():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000720725_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.231625587176232
    f_true = [[-0.27930711, -0.58474708, -0.06807932],
              [-0.01960825,  0.44659068,  0.02440941],
              [ 0.13613458, -0.38552375,  0.07408101],
              [ 0.16278078,  0.52368015, -0.0304111 ]]
    s_true = [9.69427897, 8.84439515, 8.72884187, 0.03368536, -0.22656027, 0.15691132]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval9():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000820734_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.232726625053889
    f_true = [[-0.28120461, -0.58521213, -0.06815894],
              [-0.01603377,  0.44708881,  0.02468302],
              [ 0.13491204, -0.3816136 ,  0.0734513 ],
              [ 0.16232634,  0.51973692, -0.02997538]]
    s_true = [9.24938139, 8.44580388, 8.33079108, 0.03567479, -0.23239605, 0.15246932]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval10():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020000910050_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.237556893680226
    f_true = [[-0.2743661 , -0.5711308 , -0.06601815],
              [-0.01623269,  0.43655356,  0.02358225],
              [ 0.13219026, -0.36316389,  0.07142167],
              [ 0.15840854,  0.49774113, -0.02898577]]
    s_true = [9.24980468, 8.42377645, 8.31594258, 0.04336532, -0.25156225, 0.16125368]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)


def test_eval11():
    """Test property calculations."""
    model = "Polymlp_Seko_2022_AgAu__MO_020001010051_000"
    atoms = get_structure1()
    energy, forces, stress = eval_using_kim(model, atoms)

    e_true = -11.25404325111238
    f_true = [[-0.27998978, -0.57098227, -0.06644362],
              [-0.01180036,  0.43236548,  0.02291581],
              [ 0.13282434, -0.36171296,  0.07206981],
              [ 0.1589658 ,  0.50032975, -0.028542  ]]
    s_true = [8.50236777, 8.2953701, 8.17003422, 0.02802328, -0.28749094, 0.1768198]
    np.testing.assert_allclose(energy, e_true, atol=1e-10)
    np.testing.assert_allclose(forces, f_true, atol=1e-7)
    np.testing.assert_allclose(stress, s_true, atol=1e-7)
