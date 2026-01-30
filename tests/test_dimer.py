"""An example to use kimpy to compute the energy vs. distance of an aluminum dimer."""

import numpy as np
import kimpy


# Callback routine for the KIM API to get atom neighbor lists
def get_neigh(data_in, cutoffs_in, neighbor_list_index_in, particle_number_in):
    neighbors = data_in["neighbors"][particle_number_in]
    return (neighbors, 0)


# Function to construct neighbor lists for the atoms.
# (This is a highly-inefficient N^2 algorithm that should only be used for
# small numbers of atoms.)
# For larger systems, you might want to use kimpy.neighlist.
def create_neigh(coords_in, cutoff, neigh_in):
    n = coords_in.shape[0]
    neighbors = []
    for i in range(n):
        neigh_i = []
        for j in range(n):
            if j == i:
                continue
            dist = np.linalg.norm(coords_in[i] - coords_in[j])
            if dist < cutoff:
                neigh_i.append(j)
        neigh_i = np.array(neigh_i, dtype=np.intc)
        neighbors.append(neigh_i)
    neigh_in["cutoff"] = cutoff
    neigh_in["num_particles"] = n
    neigh_in["neighbors"] = neighbors


# Initialize KIM potential specifying the units requested from the potential.
modelname = "Polymlp_Seko_2022_AgAu__MO_020000020182_000"
units_accepted, kim_model = kimpy.model.create(
    kimpy.numbering.zeroBased,
    kimpy.length_unit.A,
    kimpy.energy_unit.eV,
    kimpy.charge_unit.e,
    kimpy.temperature_unit.K,
    kimpy.time_unit.ps,
    modelname,
)

def test_dimer():
    # Define simulation arguments
    N = 2
    coords = np.zeros((N, 3), dtype=np.double)
    forces = np.zeros((N, 3), dtype=np.double)
    energy = np.array([0.0], dtype=np.double)
    num_particles = np.array([N], dtype=np.intc)
    species_code = np.zeros(num_particles, dtype=np.intc)
    particle_contributing = np.zeros(num_particles, dtype=np.intc)
    
    # Set KIM API pointers to simulation arguments
    compute_arguments = kim_model.compute_arguments_create()
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.numberOfParticles, num_particles
    )
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.particleSpeciesCodes, species_code
    )
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.particleContributing, particle_contributing
    )
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.coordinates, coords
    )
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.partialEnergy, energy
    )
    compute_arguments.set_argument_pointer(
        kimpy.compute_argument_name.partialForces, forces
    )
    
    # Setup neighbor lists through KIM API callback mechanism
    neigh = dict()
    compute_arguments.set_callback(
        kimpy.compute_callback_name.GetNeighborList, get_neigh, neigh
    )
    
    # Setup for calculation
    influence_dist = kim_model.get_influence_distance()
    supported, code = kim_model.get_species_support_and_code(kimpy.species_name.Ag)
    species_code[:] = code
    particle_contributing[:] = 1
    
    # Loop over dimer energy
    #print("  Distance", " " * 9, "Energy")
    distance_energy = []
    for z in np.arange(1.0, np.ceil(influence_dist) + 0.5, 0.5):
        coords[1, 2] = z
        create_neigh(coords, influence_dist, neigh)
        kim_model.compute(compute_arguments)
        # print("{:18.10e} {:18.10e}".format(z, energy[0]))
        distance_energy.append([z, energy[0]])

    distance_energy_true = [
        [1.0000000000e+00,  8.4606845162e+01],
        [1.5000000000e+00,  1.8009852928e+01],
        [2.0000000000e+00,  1.8196348697e+00],
        [2.5000000000e+00, -1.0521710605e+00],
        [3.0000000000e+00, -9.6529769547e-01],
        [3.5000000000e+00, -6.3299819918e-01],
        [4.0000000000e+00, -3.8103275591e-01],
        [4.5000000000e+00, -2.0442800408e-01],
        [5.0000000000e+00, -7.2369153013e-02],
        [5.5000000000e+00, -1.0103716587e-02],
        [6.0000000000e+00,  0.0000000000e+00],
    ]
    np.testing.assert_allclose(distance_energy, distance_energy_true)
