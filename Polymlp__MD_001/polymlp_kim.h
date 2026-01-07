/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#define POLYMLP_KIM

#include "polymlp_mlpcpp.h"
#include "polymlp_structs.h"
#include "polymlp_api.h"
#include "polymlp_functions_interface.h"


class PolymlpKIM {

    PolymlpAPI polymlp;
    double cutmax;
    vector1i types;

    void compute_pair(int eflag, int vflag);
    void compute_gtinv(int eflag, int vflag);

    // for pair
    void compute_antp(vector2d& antp);
    void compute_sum_of_prod_antp(
        const vector2d& antp, 
        vector2d& prod_sum_e, 
        vector2d& prod_sum_f
    );

    // for gtinv
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

    public:
    PolymlpKIM(class LAMMPS *);
    virtual ~PolymlpKIM();
    void compute(int, int);

    void parse_polymlp(int, char **);
 
};

#endif

