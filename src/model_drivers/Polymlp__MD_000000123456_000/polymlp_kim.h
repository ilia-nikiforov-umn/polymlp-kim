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

#define DIMENSION 3

typedef double VectorOfSizeDIM[DIMENSION];
typedef double VectorOfSizeSix[6];


class PolymlpKIM {

    PolymlpAPI polymlp;
    double cutoff;
    int n_atoms_contrib;

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
        VectorOfSizeSix *& particle_virial,
        bool compute_process_dEdr);

    void compute_anlmtp(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2dc& anlmtp);

    void compute_sum_of_prod_anlmtp(
        const vector2dc& anlmtp, 
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        vector2dc& prod_sum_e, 
        vector2dc& prod_sum_f
    );

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
        VectorOfSizeSix *& particle_virial,
        bool compute_process_dEdr);

    void compute_antp(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const VectorOfSizeDIM *& atom_coords,
        vector2d& antp);

    void compute_sum_of_prod_antp(
        const vector2d& antp, 
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f
    );

    public:

    PolymlpKIM(
        const std::string& parse_polymlp,
        double energy_conv,
        double, 
        double length_conv,
        double inv_length_conv,
        double);
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
        VectorOfSizeSix *& particle_virial,
        bool compute_process_dEdr
    );

    double const * cutoff_ptr() const {
        return &cutoff;
    }
 
};

#endif
