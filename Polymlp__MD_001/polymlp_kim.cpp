/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Atsuto Seko
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "polymlp_kim.h"

using namespace model_driver_Tersoff;

//#include <omp.h>

PolymlpKIM::PolymlpKIM(){}

PolymlpKIM::PolymlpKIM(
    const std::string& polymlp_file,
    // Conversion factors.
    double energy_conv,
    double, // unused inv_energy_conv
    double length_conv,
    double inv_length_conv,
    double // unused charge_conv
){
    parse_polymlp(polymlp_file, energy_conv, length_conv, inv_length_conv);

    // TODO: Is to_spec needed?
    const int n_spec = ele_strings.size();
    for (int i = 0; i < n_spec; ++i){
        to_spec[i] = ele_strings[i];
    }
}

PolymlpKIM::~PolymlpKIM(){}


void PolymlpKIM::parse_polymlp(
    const std::string& polymlp_file,
    double energy_conv,
    double length_conv,
    double inv_length_conv)
{
    // Parse polymlp file.

    polymlp.parse_polymlp_file(polymlp_file.c_str(), ele_strings, mass);

    // TODO: Unit conversion.

    const auto& fp = polymlp.get_fp();
    cutoff = fp.cutoff;

    // TODO: How to set mass values ? 

    /*
    for (int i = 1; i <= atom->ntypes; ++i){
        atom->set_mass(FLERR,i,mass[map[i-1]]);
        for (int j = 1; j <= atom->ntypes; ++j) setflag[i][j] = 1;
    }
    for (int i = 0; i < atom->natoms; ++i){
        types.emplace_back(map[(atom->type)[i]-1]);
    }
    */
}

void PolymlpKIM::compute(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms (including ghost atoms?)
    const int * const atom_types,
    const int * const contributing,
    const Array2D<const double>& atom_coords,
    double* energy, 
    double* atom_energy,
    Array2D<double>* forces,
    double* virial,
    Array2D<double>* particle_virial,
    bool compute_process_dEdr
){
    // Compute properties.
    const auto& fp = polymlp.get_fp();

    // If requested, reset energy.
    if (energy)
        *energy = 0.0;
    if (atom_energy)
        for (int i = 0; i != n_atoms; ++i) {
            atom_energy[i] = 0.0;
        }

    // Reset forces.
    if (forces)
        for (int i = 0; i != n_atoms; ++i) {
            (*forces)(i, 0) = 0.0;
            (*forces)(i, 1) = 0.0;
            (*forces)(i, 2) = 0.0;
        }

    // Reset virial.
    if (virial)
        for (int i = 0; i != 6; ++i)
            virial[i] = 0.0;
    if (particle_virial)
        for (int i = 0; i != n_atoms; ++i)
            for (int j = 0; j != 6; ++j)
                (*particle_virial)(i, j) = 0.0;

    //if (fp.feature_type == "pair"){
    //    compute_pair();
    //}

    n_atoms_contrib = 0;
    for (int i = 0; i != n_atoms; ++i) {
        if (contributing[i]) 
            ++n_atoms_contrib;
    }

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
}


void PolymlpKIM::compute_anlmtp(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms (including ghost atoms?)
    const int * const atom_types,
    const int * const contributing,
    const Array2D<const double>& atom_coords,
    vector2dc& anlmtp){

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;

    vector2d anlmtp_r(n_atoms_contrib), anlmtp_i(n_atoms_contrib);

    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.

    for (int i = 0; i != n_atoms; ++i) {
    // Skip central ghost atoms.
        if (!contributing[i]) continue;

        int itype = atom_types[i];
        const auto& maps_type = maps.maps_type[itype];
        const auto& nlmtp_attrs_noconj = maps_type.nlmtp_attrs_noconj;
        // TODO: How to specify only indices of non-ghost atoms
        anlmtp_r[i] = vector1d(nlmtp_attrs_noconj.size(), 0.0);
        anlmtp_i[i] = vector1d(nlmtp_attrs_noconj.size(), 0.0);
    }

    /*
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    */
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
        const double xtmp = atom_coords(i, 0);
        const double ytmp = atom_coords(i, 1);
        const double ztmp = atom_coords(i, 2);

        // int i,j,type1,type2,tp,jnum,*ilist,*jlist;
        int tp;
        double delx, dely, delz, dis;
        const auto& maps_type = maps.maps_type[itype];
        const auto& nlmtp_attrs_noconj = maps_type.nlmtp_attrs_noconj;

        vector1d fn; vector1dc ylm; dc val;
        for (int jj = 0; jj != n_neigh; ++jj) {
            int j = neighbors[jj];
            // TODO: Check neighbor atoms in half-neighbor list.
            // delx = x[i][0]-x[j][0];
            // dely = x[i][1]-x[j][1];
            // delz = x[i][2]-x[j][2];
            delx = xtmp - atom_coords(j,0);
            dely = ytmp - atom_coords(j,1);
            delz = ztmp - atom_coords(j,2);

            dis = sqrt(delx*delx + dely*dely + delz*delz);
            std::cout << dis << std::endl;
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
                        const int idx_j = nlmtp.jlocal_noconj_id;
                        val = fn[nlmtp.n_id] * ylm[lm_attr.ylmkey];
                        /*
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        */
                        anlmtp_r[i][idx_i] += val.real();
                        /*
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        */
                        anlmtp_r[j][idx_j] += val.real() * lm_attr.sign_j;
                        /*
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        */
                        anlmtp_i[i][idx_i] += val.imag();
                        /*
                        #ifdef _OPENMP
                        #pragma omp atomic
                        #endif
                        */
                        anlmtp_i[j][idx_j] += val.imag() * lm_attr.sign_j;
                    }
                }
            }
        }
    }
    compute_anlmtp_conjugate(
        n_atoms, atom_types, contributing, anlmtp_r, anlmtp_i, anlmtp);
}

void PolymlpKIM::compute_anlmtp_conjugate(
    int n_atoms, // Actual number of atoms (including ghost atoms?)
    const int * const atom_types,
    const int * const contributing,
    const vector2d& anlmtp_r, 
    const vector2d& anlmtp_i, 
    vector2dc& anlmtp
){

    anlmtp = vector2dc(n_atoms_contrib);

    /*
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    */
    for (int i = 0; i != n_atoms; ++i) {
        if (!contributing[i]) continue;
        const int itype = atom_types[i];
        polymlp.compute_anlmtp_conjugate(anlmtp_r[i], anlmtp_i[i], itype, anlmtp[i]);
    }
}


void PolymlpKIM::compute_gtinv(
    const KIM::ModelComputeArguments& model_compute_arguments,
    int n_atoms, // Actual number of atoms
    const int * const atom_types,
    const int * const contributing,
    const Array2D<const double>& atom_coords,
    double* energy,
    double* atom_energy,
    Array2D<double>* forces,
    double* virial,
    Array2D<double>* particle_virial,
    bool compute_process_dEdr)
{
    // Compute properties using polymlp with polynomial invariants.
    int error;       // KIM error code.
    int n_neigh;     // Number of neighbors of i.
    const int * neighbors;  // The indices of the neighbors.
    const bool eflag = energy || atom_energy; // Calculate energy?

    // loop over full neighbor list of my atoms

    vector2dc anlmtp, prod_sum_e, prod_sum_f;
    compute_anlmtp(
        model_compute_arguments,
        n_atoms,
        atom_types,
        contributing,
        atom_coords,
        anlmtp);

/*
    compute_sum_of_prod_anlmtp(anlmtp, prod_sum_e, prod_sum_f);

    vector2d evdwl_array(inum), fx_array(inum), fy_array(inum), fz_array(inum);
    for (int ii = 0; ii < inum; ii++) {
        int i = list->ilist[ii];
        int jnum = list->numneigh[i];
        evdwl_array[ii].resize(jnum);
        fx_array[ii].resize(jnum);
        fy_array[ii].resize(jnum);
        fz_array[ii].resize(jnum);
    }

    const auto& fp = polymlp.get_fp();
    const auto& maps = polymlp.get_maps();
    const auto& type_pairs = maps.type_pairs;
    const auto& tp_to_params = maps.tp_to_params;
*/
    /*
    #ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
    #endif
    */
/*
    for (int ii = 0; ii < inum; ii++) {
        int i,j,jnum,*jlist,type1,type2,tp,tagi,tagj;
        double delx,dely,delz,dis,evdwl,fx,fy,fz;
        dc val,valx,valy,valz,d1;
        vector1d fn,fn_d;
        vector1dc ylm,ylm_dx,ylm_dy,ylm_dz;

        double **x = atom->x;
        tagint *tag = atom->tag;

        i = list->ilist[ii];
        tagi = tag[i]-1;
        type1 = types[tagi];
        jnum = list->numneigh[i];
        jlist = list->firstneigh[i];

        const auto& maps_type = maps.maps_type[type1];
        const auto& nlmtp_attrs_noconj = maps_type.nlmtp_attrs_noconj;

        for (int jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            tagj = tag[j]-1;
            type2 = types[tagj];
            delx = x[i][0]-x[j][0];
            dely = x[i][1]-x[j][1];
            delz = x[i][2]-x[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                tp = type_pairs[type1][type2];
                const auto& params = tp_to_params[tp];
                const vector1d diff = {delx,dely,delz};
                const vector1d &sph = cartesian_to_spherical_(diff);
                get_fn_(dis, fp, params, fn, fn_d);
                get_ylm_(dis, sph[0], sph[1], fp.maxl, 
                         ylm, ylm_dx, ylm_dy, ylm_dz);

                evdwl = 0.0, fx = 0.0, fy = 0.0, fz = 0.0;
                for (const auto& nlmtp: nlmtp_attrs_noconj){
                    if (tp == nlmtp.tp){
                        const auto& lm_attr = nlmtp.lm;
                        const int ylmkey = lm_attr.ylmkey;
                        const int idx_i = nlmtp.ilocal_noconj_id;
                        const int idx_j = nlmtp.jlocal_noconj_id;
                        val = fn[nlmtp.n_id] * ylm[ylmkey];
                        d1 = fn_d[nlmtp.n_id] * ylm[ylmkey] / dis;
                        valx = - (d1 * delx + fn[nlmtp.n_id] * ylm_dx[ylmkey]);
                        valy = - (d1 * dely + fn[nlmtp.n_id] * ylm_dy[ylmkey]);
                        valz = - (d1 * delz + fn[nlmtp.n_id] * ylm_dz[ylmkey]);
                        const auto& prod_ei = prod_sum_e[tagi][idx_i];
                        const auto& prod_ej = prod_sum_e[tagj][idx_j];
                        const auto& prod_fi = prod_sum_f[tagi][idx_i];
                        const auto& prod_fj = prod_sum_f[tagj][idx_j];
                        const dc sum_e = prod_ei + prod_ej * lm_attr.sign_j;
                        const dc sum_f = prod_fi + prod_fj * lm_attr.sign_j;
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
                evdwl_array[ii][jj] = evdwl;
                fx_array[ii][jj] = fx;
                fy_array[ii][jj] = fy;
                fz_array[ii][jj] = fz;
            }
        }
    }

    int i,j,jnum,*jlist;
    double fx,fy,fz,evdwl,dis,delx,dely,delz;
    double **f = atom->f;
    double **x = atom->x;
    for (int ii = 0; ii < inum; ii++) {
        i = list->ilist[ii];
        jnum = list->numneigh[i], jlist = list->firstneigh[i];
        for (int jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            delx = x[i][0]-x[j][0];
            dely = x[i][1]-x[j][1];
            delz = x[i][2]-x[j][2];
            dis = sqrt(delx*delx + dely*dely + delz*delz);
            if (dis < fp.cutoff){
                evdwl = evdwl_array[ii][jj];
                fx = fx_array[ii][jj]; 
                fy = fy_array[ii][jj]; 
                fz = fz_array[ii][jj]; 
                f[i][0] += fx, f[i][1] += fy, f[i][2] += fz;
                // if (newton_pair || j < nlocal)
                f[j][0] -= fx, f[j][1] -= fy, f[j][2] -= fz;
                if (evflag) {
                    ev_tally_xyz(i,j,nlocal,newton_pair,
                            evdwl,0.0,fx,fy,fz,delx,dely,delz);
                }
            }
        }
    }
*/
}


// void PolymlpKIM::compute_sum_of_prod_anlmtp(
//     const vector2dc& anlmtp, vector2dc& prod_sum_e, vector2dc& prod_sum_f
// ){
// 
//     const int inum = list->inum;
//     prod_sum_e = vector2dc(inum);
//     prod_sum_f = vector2dc(inum);
// 
//     /*
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(guided)
//     #endif
//     */
//     for (int ii = 0; ii < inum; ii++) {
//         tagint *tag = atom->tag;
//         const int i = list->ilist[ii];
//         const int type1 = types[tag[i]-1];
//         polymlp.compute_sum_of_prod_anlmtp(
//             anlmtp[tag[i]-1], type1, prod_sum_e[tag[i]-1], prod_sum_f[tag[i]-1]
//         );
//     }
// }


// void PolymlpKIM::compute_pair(){
//     /* Compute properties using polymlp only with pair features. */
// 
//     int inum = list->inum;
//     int nlocal = atom->nlocal;
//     int newton_pair = force->newton_pair;
// 
//     vector2d antp, prod_sum_e, prod_sum_f;
//     compute_antp(antp);
//     compute_sum_of_prod_antp(antp, prod_sum_e, prod_sum_f);
// 
//     vector2d evdwl_array(inum), fpair_array(inum);
//     for (int ii = 0; ii < inum; ii++) {
//         int i = list->ilist[ii];
//         int jnum = list->numneigh[i];
//         evdwl_array[ii].resize(jnum);
//         fpair_array[ii].resize(jnum);
//     }
// 
//     const auto& fp = polymlp.get_fp();
//     const auto& maps = polymlp.get_maps();
//     const auto& type_pairs = maps.type_pairs;
//     const auto& tp_to_params = maps.tp_to_params;
// 
//     /*
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(guided)
//     #endif
//     */
//     for (int ii = 0; ii < inum; ii++) {
//         int i,j,jnum,*jlist,type1,type2,tp,tagi,tagj;
//         double delx,dely,delz,dis,evdwl,fpair;
// 
//         double **x = atom->x;
//         tagint *tag = atom->tag;
// 
//         i = list->ilist[ii];
//         tagi = tag[i]-1;
//         type1 = types[tagi];
//         jnum = list->numneigh[i];
//         jlist = list->firstneigh[i];
// 
//         const auto& maps_type = maps.maps_type[type1];
//         const auto& ntp_attrs = maps_type.ntp_attrs;
// 
//         vector1d fn,fn_d;
//         for (int jj = 0; jj < jnum; jj++) {
//             j = jlist[jj];
//             tagj = tag[j]-1;
//             type2 = types[tagj];
//             delx = x[i][0]-x[j][0];
//             dely = x[i][1]-x[j][1];
//             delz = x[i][2]-x[j][2];
//             dis = sqrt(delx*delx + dely*dely + delz*delz);
//             if (dis < fp.cutoff){
//                 tp = type_pairs[type1][type2];
//                 const auto& params = tp_to_params[tp];
//                 get_fn_(dis, fp, params, fn, fn_d);
//                 evdwl = 0.0, fpair = 0.0;
//                 for (const auto& ntp: ntp_attrs){
//                     if (tp == ntp.tp){
//                         const int idx_i = ntp.ilocal_id;
//                         const int idx_j = ntp.jlocal_id;
//                         const auto& prod_ei = prod_sum_e[tagi][idx_i];
//                         const auto& prod_ej = prod_sum_e[tagj][idx_j];
//                         const auto& prod_fi = prod_sum_f[tagi][idx_i];
//                         const auto& prod_fj = prod_sum_f[tagj][idx_j];
//                         evdwl += fn[ntp.n_id] * (prod_ei + prod_ej);
//                         fpair += fn_d[ntp.n_id] * (prod_fi + prod_fj);
//                     }
//                 }
//                 fpair *= - 1.0 / dis;
//                 evdwl_array[ii][jj] = evdwl;
//                 fpair_array[ii][jj] = fpair;
//             }
//         }
//     }
// 
//     int i,j,jnum,*jlist;
//     double fpair,evdwl,dis,delx,dely,delz;
//     double **f = atom->f;
//     double **x = atom->x;
//     for (int ii = 0; ii < inum; ii++) {
//         i = list->ilist[ii];
//         jnum = list->numneigh[i], jlist = list->firstneigh[i];
//         for (int jj = 0; jj < jnum; jj++) {
//             j = jlist[jj];
//             delx = x[i][0]-x[j][0];
//             dely = x[i][1]-x[j][1];
//             delz = x[i][2]-x[j][2];
//             dis = sqrt(delx*delx + dely*dely + delz*delz);
//             if (dis < fp.cutoff){
//                 evdwl = evdwl_array[ii][jj];
//                 fpair = fpair_array[ii][jj];
//                 f[i][0] += fpair*delx;
//                 f[i][1] += fpair*dely;
//                 f[i][2] += fpair*delz;
//                 f[j][0] -= fpair*delx;
//                 f[j][1] -= fpair*dely;
//                 f[j][2] -= fpair*delz;
//                 if (evflag) {
//                     ev_tally(i,j,nlocal,newton_pair,
//                             evdwl,0.0,fpair,delx,dely,delz);
//                 }
//             }
//         }
//     }
// }
// 
// void PolymlpKIM::compute_antp(vector2d& antp){
// 
//     const auto& fp = polymlp.get_fp();
//     const auto& maps = polymlp.get_maps();
//     const auto& type_pairs = maps.type_pairs;
//     const auto& tp_to_params = maps.tp_to_params;
// 
//     int inum = list->inum;
//     antp = vector2d(inum);
//     for (int ii = 0; ii < inum; ii++){
//         tagint *tag = atom->tag;
//         int i = list->ilist[ii];
//         int type1 = types[tag[i]-1];
// 
//         const auto& maps_type = maps.maps_type[type1];
//         const auto& ntp_attrs = maps_type.ntp_attrs;
//         antp[tag[i]-1] = vector1d(ntp_attrs.size(), 0.0);
//     }
// 
//     /*
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(auto)
//     #endif
//     */
//     for (int ii = 0; ii < inum; ii++) {
//         int i,j,type1,type2,tp,jnum,*ilist,*jlist;
//         double delx,dely,delz,dis;
// 
//         double **x = atom->x;
//         tagint *tag = atom->tag;
// 
//         i = list->ilist[ii];
//         type1 = types[tag[i]-1];
//         jnum = list->numneigh[i];
//         jlist = list->firstneigh[i];
// 
//         const auto& maps_type = maps.maps_type[type1];
//         const auto& ntp_attrs = maps_type.ntp_attrs;
// 
//         vector1d fn; 
//         for (int jj = 0; jj < jnum; ++jj) {
//             j = jlist[jj];
//             delx = x[i][0]-x[j][0];
//             dely = x[i][1]-x[j][1];
//             delz = x[i][2]-x[j][2];
//             dis = sqrt(delx*delx + dely*dely + delz*delz);
//             if (dis < fp.cutoff){
//                 type2 = types[tag[j]-1];
//                 tp = type_pairs[type1][type2];
//                 const auto& params = tp_to_params[tp];
//                 get_fn_(dis, fp, params, fn);
//                 for (const auto& ntp: ntp_attrs){
//                     if (tp == ntp.tp){
//                         const int idx_i = ntp.ilocal_id;
//                         const int idx_j = ntp.jlocal_id;
//                         /*
//                         #ifdef _OPENMP
//                         #pragma omp atomic
//                         #endif
//                         */
//                         antp[tag[i]-1][idx_i] += fn[ntp.n_id];
//                         /*
//                         #ifdef _OPENMP
//                         #pragma omp atomic
//                         #endif
//                         */
//                         antp[tag[j]-1][idx_j] += fn[ntp.n_id];
//                     }
//                 }
//             }
//         }
//     }
// }
// 
// void PolymlpKIM::compute_sum_of_prod_antp(
//     const vector2d& antp, vector2d& prod_sum_e, vector2d& prod_sum_f
// ){
//     const int inum = list->inum;
//     prod_sum_e = vector2d(inum);
//     prod_sum_f = vector2d(inum);
// 
//     /*
//     #ifdef _OPENMP
//     #pragma omp parallel for schedule(guided)
//     #endif
//     */
//     for (int ii = 0; ii < inum; ii++) {
//         tagint *tag = atom->tag;
//         const int i = list->ilist[ii];
//         const int type1 = types[tag[i]-1];
//         polymlp.compute_sum_of_prod_antp(
//             antp[tag[i]-1], type1, prod_sum_e[tag[i]-1], prod_sum_f[tag[i]-1]
//         );
//     }
// }

/* ---------------------------------------------------------------------- */


