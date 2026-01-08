/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#define POLYMLP_KIM

#include "polymlp/polymlp_mlpcpp.h"
#include "polymlp/polymlp_structs.h"
#include "polymlp/polymlp_api.h"
#include "polymlp/polymlp_functions_interface.h"

#include "KIM_ModelDriverHeaders.hpp"
#include "ndarray.hpp"


class PolymlpKIM {

    PolymlpAPI polymlp;
    double cutoff;
    vector1i types;

    int n_atoms_contrib;

    // Parse polymlp file.
    std::vector<std::string> ele_strings;
    vector1d mass;
    void parse_polymlp(
        const std::string& parse_polymlp,
        double energy_conv,
        double length_conv,
        double inv_length_conv);

    // Compute properties using polymlp with pairwise features.
    void compute_pair();
    void compute_antp(vector2d& antp);
    void compute_sum_of_prod_antp(
        const vector2d& antp, 
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f
    );

    // Compute properties using polymlp with polynomial invariants.
    void compute_gtinv();
    void compute_anlmtp(vector2dc& anlmtp);
    void compute_anlmtp_conjugate(
        const vector2d& anlmtp_r, 
        const vector2d& anlmtp_i, 
        vector2dc& anlmtp
    );
    void compute_sum_of_prod_anlmtp(
        const vector2dc& anlmtp, 
        vector2dc& prod_sum_e, 
        vector2dc& prod_sum_f
    );

    protected:

    // map element index to element name, needed for user-friendly error messages
    std::map<int, std::string> to_spec;

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
    void compute();
 
};

#endif

