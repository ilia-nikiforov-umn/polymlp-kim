/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#ifndef POLYMLP_KIM
#define POLYMLP_KIM

#include "polymlp/polymlp_mlpcpp.h"
#include "polymlp/polymlp_structs.h"
#include "polymlp/polymlp_api.h"
#include "polymlp/polymlp_functions_interface.h"

#include "KIM_ModelDriverHeaders.hpp"

#include <omp.h>

#define DIMENSION 3

typedef double VectorOfSizeDIM[DIMENSION];
typedef double VectorOfSizeSix[6];


class PolymlpKIM {

    std::vector<PolymlpAPI> polymlp_array;
    PolymlpAPI polymlp;
    double cutoff_max;
    int n_atoms;
    int n_atoms_contrib;

    vector1i map_contrib_to_full;
    std::map<int, int> map_full_to_contrib;

    // Compute properties using polymlp with polynomial invariants.
    void compute_gtinv(
        const KIM::ModelComputeArguments& model_compute_arguments,
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        double* energy,
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);

    void compute_sum_of_prod_anlmtp(
        const KIM::ModelComputeArguments& model_compute_arguments,
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2dc& prod_sum_e, 
        vector2dc& prod_sum_f);

    void accumulate_properties(
        const vector1d& energy_array,
        const vector2d& atom_energy_array,
        const vector2d& fx_array,
        const vector2d& fy_array,
        const vector2d& fz_array,
        const vector2d& virial_array,
        const vector3d& particle_virial_array,
        double* energy,
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);

    // Compute properties using polymlp with pairwise features.
    void compute_pair(
        const KIM::ModelComputeArguments& model_compute_arguments,
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        double* energy,
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);

    void compute_sum_of_prod_antp(
        const KIM::ModelComputeArguments& model_compute_arguments,
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f);

    public:

    PolymlpKIM(
        const std::vector<std::string>& polymlp_files,
        double energy_conv,
        double length_conv,
        double inv_length_conv);
    ~PolymlpKIM();

    void compute(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms_in, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        double* energy, 
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);

    double const * cutoff_ptr() const {
        return &cutoff_max;
    }

    int array_idx(int tid, int i) {
        return tid * n_atoms + i;
    }
};

#endif
