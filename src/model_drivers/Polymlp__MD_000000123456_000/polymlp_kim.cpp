/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#include "polymlp_kim.h"


PolymlpKIM::PolymlpKIM(
    const std::string& polymlp_file,
    // Conversion factors.
    double energy_conv,
    double, // unused inv_energy_conv
    double length_conv,
    double inv_length_conv,
    double // unused charge_conv
){

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
    VectorOfSizeSix *& particle_virial,
    bool compute_process_dEdr)
{
    // Compute properties.

    // Set number of contributing atoms.
    n_atoms_contrib = 0;
    for (int i = 0; i != n_atoms; ++i) {
        if (contributing[i]) 
            ++n_atoms_contrib;
    }
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
            particle_virial,
            compute_process_dEdr);
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
            particle_virial,
            compute_process_dEdr);
    }
}


void PolymlpKIM::compute_anlmtp(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    vector2dc& anlmtp){

    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    anlmtp = vector2dc(n_atoms_contrib);
    for (int i = 0; i != n_atoms; ++i) {
    // Skip central ghost atoms.
        if (!contributing[i]) continue;
        // Get neighbors.
        error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
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
        polymlp.compute_anlmtp_conjugate(
            anlmtp_r, anlmtp_i, atom_types[i], anlmtp[i]
        );
    }
}


void PolymlpKIM::compute_sum_of_prod_anlmtp(
    const vector2dc& anlmtp, 
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    vector2dc& prod_sum_e, 
    vector2dc& prod_sum_f){

    prod_sum_e = vector2dc(n_atoms_contrib);
    prod_sum_f = vector2dc(n_atoms_contrib);
    for (int i = 0; i != n_atoms; ++i) {
        // Skip central ghost atoms.
        if (!contributing[i]) continue;
 
        const int itype = atom_types[i];
        polymlp.compute_sum_of_prod_anlmtp(
            anlmtp[i], itype, prod_sum_e[i], prod_sum_f[i]);
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
    VectorOfSizeSix *& particle_virial,
    bool compute_process_dEdr)
{
    // Compute properties using polymlp with polynomial invariants.
    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    vector2dc anlmtp, prod_sum_e, prod_sum_f;
    compute_anlmtp(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        anlmtp);

    compute_sum_of_prod_anlmtp(
        anlmtp, 
        n_atoms,
        atom_types,
        contributing,
        prod_sum_e, 
        prod_sum_f);

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    for (int i = 0; i != n_atoms; ++i) {
        // Skip central ghost atoms.
        if (!contributing[i]) continue;

        // Get neighbors.
        error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
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
                        const auto& prod_ei = prod_sum_e[i][idx_i];
                        const auto& prod_fi = prod_sum_f[i][idx_i];
                        if (contributing[j]){
                            const int idx_j = nlmtp.jlocal_noconj_id;
                            const auto& prod_ej = prod_sum_e[j][idx_j];
                            const auto& prod_fj = prod_sum_f[j][idx_j];
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


void PolymlpKIM::compute_antp(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    const VectorOfSizeDIM *& atom_coords,
    vector2d& antp){

    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    antp = vector2d(n_atoms_contrib);
    for (int i = 0; i != n_atoms; ++i) {
    // Skip central ghost atoms.
        if (!contributing[i]) continue;
        // Get neighbors.
        error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
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

        antp[i] = vector1d(ntp_attrs.size(), 0.0);
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
                        antp[i][idx_i] += fn[ntp.n_id];
                    }
                }
            }
        }
    }
}


void PolymlpKIM::compute_sum_of_prod_antp(
    const vector2d& antp, 
    int n_atoms, // Actual number of atoms including ghost atoms
    const int * const atom_types,
    const int * const contributing,
    vector2d& prod_sum_e, vector2d& prod_sum_f
){
    prod_sum_e = vector2d(n_atoms_contrib);
    prod_sum_f = vector2d(n_atoms_contrib);

    for (int i = 0; i != n_atoms; ++i) {
        if (!contributing[i]) continue;
 
        const int itype = atom_types[i];
        polymlp.compute_sum_of_prod_antp(
            antp[i], itype, prod_sum_e[i], prod_sum_f[i]);
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
    VectorOfSizeSix *& particle_virial,
    bool compute_process_dEdr)
{
    // Compute properties using polymlp only with pair features.

    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    vector2d antp, prod_sum_e, prod_sum_f;
    compute_antp(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        antp);

    compute_sum_of_prod_antp(
        antp, 
        n_atoms,
        atom_types,
        contributing,
        prod_sum_e, 
        prod_sum_f);

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    for (int i = 0; i != n_atoms; ++i) {
        // Skip central ghost atoms.
        if (!contributing[i]) continue;

        // Get neighbors.
        error = model_compute_arguments.GetNeighborList(0, i, &n_neigh, &neighbors);
        if (error) {
            throw std::runtime_error(
                "Error in KIM::ModelComputeArguments.GetNeighborList");
        }
        const int itype = atom_types[i];
        const double xtmp = atom_coords[i][0];
        const double ytmp = atom_coords[i][1];
        const double ztmp = atom_coords[i][2];

        int tp;
        double delx,dely,delz,dis,evdwl,fpair,fx,fy,fz;
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
                        const auto& prod_ei = prod_sum_e[i][idx_i];
                        const auto& prod_fi = prod_sum_f[i][idx_i];
                        double val_e, val_f;
                        if (contributing[j]){
                            const int idx_j = ntp.jlocal_id;
                            const auto& prod_ej = prod_sum_e[j][idx_j];
                            const auto& prod_fj = prod_sum_f[j][idx_j];
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
                fx = fpair * delx;
                fy = fpair * dely;
                fz = fpair * delz;

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
