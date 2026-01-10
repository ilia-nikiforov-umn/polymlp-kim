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
#include "ndarray.hpp"


using namespace model_driver_Tersoff;

class PolymlpKIM {

    PolymlpAPI polymlp;
    double cutoff;
    int n_atoms_contrib;

    // Parse polymlp file.
    std::vector<std::string> ele_strings;
    vector1d mass;
    void parse_polymlp(
        const std::string& parse_polymlp,
        double energy_conv,
        double length_conv,
        double inv_length_conv);

    // Compute properties using polymlp with polynomial invariants.
    void compute_gtinv(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms,
        const int * const atom_types,
        const int * const contributing,
        const Array2D<const double>& atom_coords,
        double* energy,
        double* atom_energy,
        Array2D<double>* forces,
        double* virial,
        Array2D<double>* particle_virial,
        bool compute_process_dEdr);

    void compute_anlmtp(
        const KIM::ModelComputeArguments& model_compute_arguments,
        int n_atoms, 
        const int * const atom_types,
        const int * const contributing,
        const Array2D<const double>& atom_coords,
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
    /*
    void compute_pair();
    void compute_antp(vector2d& antp);
    void compute_sum_of_prod_antp(
        const vector2d& antp, 
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f
    );
    */

    // protected:

    public:

    PolymlpKIM();
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
        const Array2D<const double>& atom_coords,
        double* energy, 
        double* atom_energy,
        Array2D<double>* forces,
        double* virial,
        Array2D<double>* particle_virial,
        bool compute_process_dEdr
    );

    double const * cutoff_ptr() const {
        return &cutoff;
    }
 
};

#endif
