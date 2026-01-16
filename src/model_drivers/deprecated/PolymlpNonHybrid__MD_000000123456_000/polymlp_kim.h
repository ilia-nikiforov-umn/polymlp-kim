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

    PolymlpAPI polymlp;
    double cutoff;
    int n_atoms_contrib;
    vector1i map_contrib_to_full;
    std::map<int, int> map_full_to_contrib;

    // Compute properties using polymlp with polynomial invariants.
    void compute_gtinv(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms,
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
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2dc& prod_sum_e, 
        vector2dc& prod_sum_f);

    void accumulate_properties(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        const vector2d& energy_array,
        const vector2d& fx_array,
        const vector2d& fy_array,
        const vector2d& fz_array,
        double* energy,
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);
 
    // Compute properties using polymlp with pairwise features.
    void compute_pair(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms,
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
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f);

    public:

    PolymlpKIM(
        const std::string& parse_polymlp,
        double energy_conv,
        double length_conv,
        double inv_length_conv);
    ~PolymlpKIM();

    void compute(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        double* energy, 
        double* atom_energy,
        VectorOfSizeDIM *& forces,
        double* virial,
        VectorOfSizeSix *& particle_virial);

    double const * cutoff_ptr() const {
        return &cutoff;
    }
 
};

#endif
