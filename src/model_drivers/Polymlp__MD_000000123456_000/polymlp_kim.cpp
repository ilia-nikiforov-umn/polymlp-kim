/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#include "polymlp_kim.h"


PolymlpKIM::PolymlpKIM(
    const std::string& polymlp_file,
    // Conversion factors.
    double energy_conv,
    double length_conv,
    double inv_length_conv)
{

    // Parse polymlp file.
    std::vector<std::string> ele_strings;
    vector1d mass;
    polymlp.parse_polymlp_file(polymlp_file.c_str(), ele_strings, mass);

    // Unit conversion.
    polymlp.convert_unit(energy_conv, length_conv, inv_length_conv);

    const auto& fp = polymlp.get_fp();
    cutoff = fp.cutoff;
}


PolymlpKIM::~PolymlpKIM(){}


void PolymlpKIM::compute(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    double* energy, 
    double* atom_energy,
    VectorOfSizeDIM *& forces,
    double* virial,
    VectorOfSizeSix *& particle_virial)
{
    // Compute properties.

    // If requested, reset energy.
    if (energy){
        *energy = 0.0;
    }
    // Reset atomic energy.
    if (atom_energy){
        for (int i = 0; i != n_atoms; ++i) {
            atom_energy[i] = 0.0;
        }
    }
    // Reset forces.
    if (forces){
        for (int i = 0; i != n_atoms; ++i) {
            forces[i][0] = 0.0;
            forces[i][1] = 0.0;
            forces[i][2] = 0.0;
        }
    }
    // Reset virial.
    if (virial){
        for (int i = 0; i != 6; ++i)
            virial[i] = 0.0;
    }
    // Reset partial virial.
    if (particle_virial){
        for (int i = 0; i != n_atoms; ++i)
            for (int j = 0; j != 6; ++j)
                particle_virial[i][j] = 0.0;
    }

    // Set mapping from contributing atom IDs to full atom IDs.
    //  and mapping from full atom IDs to contributing atom IDs.
    map_contrib_to_full = vector1i({});
    map_full_to_contrib = std::map<int, int>();
    int icontrib(0);
    for (int i = 0; i != n_atoms; ++i) {
        if (contributing[i]) {
            map_contrib_to_full.emplace_back(i);
            map_full_to_contrib[i] = icontrib;
            ++icontrib;
        }
    }
    n_atoms_contrib = map_contrib_to_full.size();

    const auto& fp = polymlp.get_fp();
    if (fp.feature_type == "gtinv"){
        compute_gtinv(
            model_compute_arguments,
            n_atoms,
            atom_types,
            contributing,
            atom_coords,
            energy,
            atom_energy,
            forces,
            virial,
            particle_virial);
    }
    else if (fp.feature_type == "pair"){
        compute_pair(
            model_compute_arguments,
            n_atoms,
            atom_types,
            contributing,
            atom_coords,
            energy,
            atom_energy,
            forces,
            virial,
            particle_virial);
    }
}


void PolymlpKIM::compute_sum_of_prod_anlmtp(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    vector2dc& prod_sum_e, 
    vector2dc& prod_sum_f){

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    prod_sum_e = vector2dc(n_atoms_contrib);
    prod_sum_f = vector2dc(n_atoms_contrib);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    for (int icontrib = 0; icontrib != n_atoms_contrib; ++icontrib) {
        const int i = map_contrib_to_full[icontrib];
        // Get neighbors.
        int n_neigh;            // Number of neighbors of i.
        const int * neighbors;  // The indices of the neighbors.
        int error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
        const int itype = atom_types[i];
        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];

        int tp;
        double delx, dely, delz, dis;
        vector1d fn; vector1dc ylm; dc val;

        const auto& maps_type = maps.maps_type[itype];
        const auto& nlmtp_attrs_noconj = maps_type.nlmtp_attrs_noconj;
        vector1d anlmtp_r(nlmtp_attrs_noconj.size(), 0.0);
        vector1d anlmtp_i(nlmtp_attrs_noconj.size(), 0.0);
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            delx = xtmp - atom_coords[j][0];
            dely = ytmp - atom_coords[j][1];
            delz = ztmp - atom_coords[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                const int jtype = atom_types[j];
                tp = type_pairs[itype][jtype];
                const vector1d &sph = cartesian_to_spherical_(vector1d{delx,dely,delz});
                const auto& params = tp_to_params[tp];
                get_fn_(dis, fp, params, fn);
                get_ylm_(sph[0], sph[1], fp.maxl, ylm);
                for (const auto& nlmtp: nlmtp_attrs_noconj){
                    if (tp == nlmtp.tp){
                        const auto& lm_attr = nlmtp.lm;
                        const int idx_i = nlmtp.ilocal_noconj_id;
                        val = fn[nlmtp.n_id] * ylm[lm_attr.ylmkey];
                        anlmtp_r[idx_i] += val.real();
                        anlmtp_i[idx_i] += val.imag();
                    }
                }
            }
        }
        vector1dc anlmtp;
        polymlp.compute_anlmtp_conjugate(anlmtp_r, anlmtp_i, itype, anlmtp);
        polymlp.compute_sum_of_prod_anlmtp(
            anlmtp, itype, prod_sum_e[icontrib], prod_sum_f[icontrib]);
    }
}


void PolymlpKIM::compute_gtinv(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    double* energy,
    double* atom_energy,
    VectorOfSizeDIM *& forces,
    double* virial,
    VectorOfSizeSix *& particle_virial)
{
    // Compute properties using polymlp with polynomial invariants.
    vector2dc anlmtp, prod_sum_e, prod_sum_f;
    compute_sum_of_prod_anlmtp(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        prod_sum_e,
        prod_sum_f);

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    // TODO: More efficient implementation using thread local arrays.
    vector2d energy_array(n_atoms_contrib);
    vector2d fx_array(n_atoms_contrib);
    vector2d fy_array(n_atoms_contrib);
    vector2d fz_array(n_atoms_contrib);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    for (int icontrib = 0; icontrib != n_atoms_contrib; ++icontrib) {
        const int i = map_contrib_to_full[icontrib];
        // Get neighbors.
        int n_neigh;            // Number of neighbors of i.
        const int * neighbors;  // The indices of the neighbors.
        int error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
        energy_array[icontrib].resize(n_neigh);
        fx_array[icontrib].resize(n_neigh);
        fy_array[icontrib].resize(n_neigh);
        fz_array[icontrib].resize(n_neigh);

        const int itype = atom_types[i];
        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];

        int tp;
        double delx,dely,delz,dis,evdwl,fx,fy,fz;
        dc val,valx,valy,valz,d1;
        vector1d fn,fn_d;
        vector1dc ylm,ylm_dx,ylm_dy,ylm_dz;

        const auto& maps_type = maps.maps_type[itype];
        const auto& nlmtp_attrs_noconj = maps_type.nlmtp_attrs_noconj;
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            delx = xtmp - atom_coords[j][0];
            dely = ytmp - atom_coords[j][1];
            delz = ztmp - atom_coords[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                const int jtype = atom_types[j];
                tp = type_pairs[itype][jtype];
                const auto& params = tp_to_params[tp];
                const vector1d diff = {delx,dely,delz};
                const vector1d &sph = cartesian_to_spherical_(diff);
                get_fn_(dis, fp, params, fn, fn_d);
                get_ylm_(dis, sph[0], sph[1], fp.maxl, ylm, ylm_dx, ylm_dy, ylm_dz);

                evdwl = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;
                for (const auto& nlmtp: nlmtp_attrs_noconj){
                    if (tp == nlmtp.tp){
                        const auto& lm_attr = nlmtp.lm;
                        const int ylmkey = lm_attr.ylmkey;
                        val = fn[nlmtp.n_id] * ylm[ylmkey];
                        d1 = fn_d[nlmtp.n_id] * ylm[ylmkey] / dis;
                        valx = - (d1 * delx + fn[nlmtp.n_id] * ylm_dx[ylmkey]);
                        valy = - (d1 * dely + fn[nlmtp.n_id] * ylm_dy[ylmkey]);
                        valz = - (d1 * delz + fn[nlmtp.n_id] * ylm_dz[ylmkey]);

                        dc sum_e, sum_f;
                        const int idx_i = nlmtp.ilocal_noconj_id;
                        const auto& prod_ei = prod_sum_e[icontrib][idx_i];
                        const auto& prod_fi = prod_sum_f[icontrib][idx_i];
                        if (contributing[j]){
                            const int jcontrib = map_full_to_contrib[j];
                            const int idx_j = nlmtp.jlocal_noconj_id;
                            const auto& prod_ej = prod_sum_e[jcontrib][idx_j];
                            const auto& prod_fj = prod_sum_f[jcontrib][idx_j];
                            sum_e = 0.5 * (prod_ei + prod_ej * lm_attr.sign_j);
                            sum_f = 0.5 * (prod_fi + prod_fj * lm_attr.sign_j);
                        }
                        else {
                            sum_e = prod_ei;
                            sum_f = prod_fi;
                        }
                        if (lm_attr.m == 0){
                            evdwl += 0.5 * prod_real(val, sum_e);
                            fx += 0.5 * prod_real(valx, sum_f);
                            fy += 0.5 * prod_real(valy, sum_f);
                            fz += 0.5 * prod_real(valz, sum_f);
                        }
                        else {
                            evdwl += prod_real(val, sum_e);
                            fx += prod_real(valx, sum_f);
                            fy += prod_real(valy, sum_f);
                            fz += prod_real(valz, sum_f);
                        }
                    }
                }
                energy_array[icontrib][jj] = evdwl;
                fx_array[icontrib][jj] = fx;
                fy_array[icontrib][jj] = fy;
                fz_array[icontrib][jj] = fz;
            }
        }
    }

    accumulate_properties(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        energy_array,
        fx_array,
        fy_array,
        fz_array,
        energy,
        atom_energy,
        forces,
        virial,
        particle_virial);
}


void PolymlpKIM::compute_sum_of_prod_antp(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    vector2d& prod_sum_e, 
    vector2d& prod_sum_f){

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    prod_sum_e = vector2d(n_atoms_contrib);
    prod_sum_f = vector2d(n_atoms_contrib);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    for (int icontrib = 0; icontrib != n_atoms_contrib; ++icontrib) {
        const int i = map_contrib_to_full[icontrib];
        // Get neighbors.
        int n_neigh;            // Number of neighbors of i.
        const int * neighbors;  // The indices of the neighbors.
        int error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
        const int itype = atom_types[i];
        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];

        int tp;
        double delx, dely, delz, dis;
        vector1d fn;

        const auto& maps_type = maps.maps_type[itype];
        const auto& ntp_attrs = maps_type.ntp_attrs;

        vector1d antp(ntp_attrs.size(), 0.0);
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            delx = xtmp - atom_coords[j][0];
            dely = ytmp - atom_coords[j][1];
            delz = ztmp - atom_coords[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                const int jtype = atom_types[j];
                tp = type_pairs[itype][jtype];
                const auto& params = tp_to_params[tp];
                get_fn_(dis, fp, params, fn);
                for (const auto& ntp: ntp_attrs){
                    if (tp == ntp.tp){
                        const int idx_i = ntp.ilocal_id;
                        antp[idx_i] += fn[ntp.n_id];
                    }
                }
            }
        }
        polymlp.compute_sum_of_prod_antp(
            antp, itype, prod_sum_e[icontrib], prod_sum_f[icontrib]);
    }
}

void PolymlpKIM::compute_pair(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    double* energy,
    double* atom_energy,
    VectorOfSizeDIM *& forces,
    double* virial,
    VectorOfSizeSix *& particle_virial)
{
    // Compute properties using polymlp only with pair features.

    vector2d prod_sum_e, prod_sum_f;
    compute_sum_of_prod_antp(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        prod_sum_e, 
        prod_sum_f);

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    vector2d energy_array(n_atoms_contrib);
    vector2d fx_array(n_atoms_contrib);
    vector2d fy_array(n_atoms_contrib);
    vector2d fz_array(n_atoms_contrib);
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    for (int icontrib = 0; icontrib != n_atoms_contrib; ++icontrib) {
        const int i = map_contrib_to_full[icontrib];
        // Get neighbors.
        int n_neigh;            // Number of neighbors of i.
        const int * neighbors;  // The indices of the neighbors.
        int error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
        energy_array[icontrib].resize(n_neigh);
        fx_array[icontrib].resize(n_neigh);
        fy_array[icontrib].resize(n_neigh);
        fz_array[icontrib].resize(n_neigh);

        const int itype = atom_types[i];
        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];

        int tp;
        double delx,dely,delz,dis,evdwl,fpair;
        vector1d fn,fn_d;

        const auto& maps_type = maps.maps_type[itype];
        const auto& ntp_attrs = maps_type.ntp_attrs;
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            delx = xtmp - atom_coords[j][0];
            dely = ytmp - atom_coords[j][1];
            delz = ztmp - atom_coords[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                const int jtype = atom_types[j];
                tp = type_pairs[itype][jtype];
                const auto& params = tp_to_params[tp];
                get_fn_(dis, fp, params, fn, fn_d);
                evdwl = 0.0, fpair = 0.0;
                for (const auto& ntp: ntp_attrs){
                    if (tp == ntp.tp){
                        const int idx_i = ntp.ilocal_id;
                        const auto& prod_ei = prod_sum_e[icontrib][idx_i];
                        const auto& prod_fi = prod_sum_f[icontrib][idx_i];
                        double val_e, val_f;
                        if (contributing[j]){
                            const int jcontrib = map_full_to_contrib[j];
                            const int idx_j = ntp.jlocal_id;
                            const auto& prod_ej = prod_sum_e[jcontrib][idx_j];
                            const auto& prod_fj = prod_sum_f[jcontrib][idx_j];
                            val_e = 0.5 * (prod_ei + prod_ej);
                            val_f = 0.5 * (prod_fi + prod_fj);
                        }
                        else {
                            val_e = prod_ei;
                            val_f = prod_fi;
                        }
                        evdwl += fn[ntp.n_id] * val_e;
                        fpair += fn_d[ntp.n_id] * val_f;
                    }
                }
                fpair *= - 1.0 / dis;
                energy_array[icontrib][jj] = evdwl;
                fx_array[icontrib][jj] = fpair * delx;
                fy_array[icontrib][jj] = fpair * dely;
                fz_array[icontrib][jj] = fpair * delz;
            }
        }
    }
    accumulate_properties(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        energy_array,
        fx_array,
        fy_array,
        fz_array,
        energy,
        atom_energy,
        forces,
        virial,
        particle_virial);

}


void PolymlpKIM::accumulate_properties(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
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
    VectorOfSizeSix *& particle_virial)
{
    int error;              // KIM error code.
    int n_neigh;            // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    for (int icontrib = 0; icontrib != n_atoms_contrib; ++icontrib) {
        const int i = map_contrib_to_full[icontrib];
        // Get neighbors.
        error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }

        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            const double delx = xtmp - atom_coords[j][0];
            const double dely = ytmp - atom_coords[j][1];
            const double delz = ztmp - atom_coords[j][2];
            const double dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < cutoff){
                const double evdwl = energy_array[icontrib][jj];
                const double fx = fx_array[icontrib][jj];
                const double fy = fy_array[icontrib][jj];
                const double fz = fz_array[icontrib][jj];
                if (energy)
                    *energy += evdwl;
                if (atom_energy){
                    atom_energy[i] += 0.5 * evdwl;
                    atom_energy[j] += 0.5 * evdwl;
                }
                if (forces){
                    forces[i][0] += fx; 
                    forces[i][1] += fy; 
                    forces[i][2] += fz;
                    forces[j][0] -= fx; 
                    forces[j][1] -= fy; 
                    forces[j][2] -= fz;
                }
                vector1d val_tmp(6);
                if (virial || particle_virial){

                    val_tmp[0] = delx * fx;
                    val_tmp[1] = dely * fy;
                    val_tmp[2] = delz * fz;
                    val_tmp[3] = dely * fz;
                    val_tmp[4] = delx * fz;
                    val_tmp[5] = delx * fy;
                    // lammps convension
                    // virial[3] += delx * fy;
                    // virial[5] += dely * fz;
                }
                if (virial){
                    for (int k = 0; k < 6; ++k){
                        virial[k] += val_tmp[k];
                    }
                }
                if (particle_virial){
                    for (int k = 0; k < 6; ++k){
                        particle_virial[i][k] += 0.5 * val_tmp[k];
                        particle_virial[j][k] += 0.5 * val_tmp[k];
                    }
                }
            }
        }
    }
}

