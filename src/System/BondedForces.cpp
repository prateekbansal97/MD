//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
#include "System/Metrics.h"

namespace md
{
    void System::calculate_forces_bonds() {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& bond: topology_.get_harmonic_bonds())
        {
            const int atomAIndex = bond.get_atomA_index();
            const int atomBIndex = bond.get_atomB_index();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            const double dx1 = x1 - x2;
            const double dy1 = y1 - y2;
            const double dz1 = z1 - z2;

            const double r = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

            const double k = bond.get_Bond_force_constant();
            const double r0 = bond.get_Bond_equil_length();

            if (r < 1e-12) continue;

            const double famag = -2 * k * (r - r0);
            const double fax = famag*dx1/r;
            const double fay = famag*dy1/r;
            const double faz = famag*dz1/r;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] -= fax; forces_[3*atomBIndex+1] -= fay; forces_[3*atomBIndex+2] -= faz;
        }
    }

    void System::calculate_forces_UBbonds()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& bond: topology_.get_harmonic_UBs())
        {
            const int atomAIndex = bond.get_atomA_index();
            const int atomBIndex = bond.get_atomB_index();
            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];

            const double dx1 = x1 - x2;
            const double dy1 = y1 - y2;
            const double dz1 = z1 - z2;

            const double r = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

            const double k = bond.get_UB_force_constant();
            const double r0 = bond.get_UB_equil_value();

            if (r < 1e-12) continue;

            const double famag = -2 * k * (r - r0);
            const double fax = famag*dx1/r;
            const double fay = famag*dy1/r;
            const double faz = famag*dz1/r;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] += -fax; forces_[3*atomBIndex+1] += -fay; forces_[3*atomBIndex+2] += -faz;
        }
    }

    void System::calculate_forces_angles()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& ang: topology_.get_harmonic_angles())
        {
            const int atomAIndex = ang.get_atomA_index();
            const int atomBIndex = ang.get_atomB_index();
            const int atomCIndex = ang.get_atomC_index();

            const double x1 = coordinates[3*atomAIndex], y1 = coordinates[3*atomAIndex + 1], z1 = coordinates[3*atomAIndex + 2];
            const double x2 = coordinates[3*atomBIndex], y2 = coordinates[3*atomBIndex + 1], z2 = coordinates[3*atomBIndex + 2];
            const double x3 = coordinates[3*atomCIndex], y3 = coordinates[3*atomCIndex + 1], z3 = coordinates[3*atomCIndex + 2];

            // ba
            const double ba_x = x1 - x2;
            const double ba_y = y1 - y2;
            const double ba_z = z1 - z2;

            // bc
            const double bc_x = x3 - x2;
            const double bc_y = y3 - y2;
            const double bc_z = z3 - z2;

            const double mod_ab_sq = ba_x*ba_x + ba_y*ba_y + ba_z*ba_z;
            const double mod_ab = std::sqrt(mod_ab_sq);

            if (mod_ab < 1e-12) continue;

            // pa = norm(ba x (ba x bc))
            //    = norm ((ba.bc)ba - (modba*modba)bc)

            const double dot_ba_bc = ba_x*bc_x + ba_y*bc_y + ba_z*bc_z;
            double pa_x = dot_ba_bc*ba_x - mod_ab_sq*bc_x;
            double pa_y = dot_ba_bc*ba_y - mod_ab_sq*bc_y;
            double pa_z = dot_ba_bc*ba_z - mod_ab_sq*bc_z;

            const double mod_pa = std::sqrt(pa_x*pa_x + pa_y*pa_y + pa_z*pa_z);
            if (mod_pa < 1e-12) continue;
            pa_x /= mod_pa;
            pa_y /= mod_pa;
            pa_z /= mod_pa;

            const double mod_bc_sq = bc_x*bc_x + bc_y*bc_y + bc_z*bc_z;
            const double mod_bc = std::sqrt(mod_bc_sq);

            if (mod_bc < 1e-12) continue;

            // pc = norm(cb x (ba x bc))
            //    = norm ((bc.ba)bc -(modbc*modbc)ba)

            double pc_x = dot_ba_bc*bc_x - mod_bc_sq*ba_x;
            double pc_y = dot_ba_bc*bc_y - mod_bc_sq*ba_y;
            double pc_z = dot_ba_bc*bc_z - mod_bc_sq*ba_z;

            const double mod_pc = std::sqrt(pc_x*pc_x + pc_y*pc_y + pc_z*pc_z);
            if (mod_pc < 1e-12) continue;
            pc_x /= mod_pc;
            pc_y /= mod_pc;
            pc_z /= mod_pc;


            const double k = ang.get_Angle_force_constant();
            const double theta0 = ang.get_Angle_equil_angle();
            const double theta = Metrics::angle(x1, y1, z1, x2, y2, z2, x3, y3, z3);

            const double fmag = -2 * k * (theta - theta0);
            const double fax = fmag/mod_ab*pa_x;
            const double fay = fmag/mod_ab*pa_y;
            const double faz = fmag/mod_ab*pa_z;

            const double fcx = fmag/mod_bc*pc_x;
            const double fcy = fmag/mod_bc*pc_y;
            const double fcz = fmag/mod_bc*pc_z;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
            forces_[3*atomBIndex] += -fax-fcx; forces_[3*atomBIndex+1] += -fay-fcy; forces_[3*atomBIndex+2] += -faz-fcz;

        }
    }

    void System::calculate_forces_cosinedihedrals() {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& dih: topology_.get_cosine_dihedrals()) {
            const int atomAIndex = dih.get_atomA_index();
            const int atomBIndex = dih.get_atomB_index();
            const int atomCIndex = dih.get_atomC_index();
            const int atomDIndex = dih.get_atomD_index();

            const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                    3 * atomAIndex + 2];
            const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                    3 * atomBIndex + 2];
            const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                    3 * atomCIndex + 2];
            const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                    3 * atomDIndex + 2];

            // ba
            const double ba_x = x1 - x2;
            const double ba_y = y1 - y2;
            const double ba_z = z1 - z2;

            // cb
            const double cb_x = x2 - x3;
            const double cb_y = y2 - y3;
            const double cb_z = z2 - z3;

            // dc
            const double dc_x = x3 - x4;
            const double dc_y = y3 - y4;
            const double dc_z = z3 - z4;


            const double k = dih.get_Dihedral_force_constant();
            const double n = dih.get_Dihedral_Periodicity();
            const double delta = dih.get_Dihedral_phase();
            const double phi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
            const double tau = k * n * std::sin(n * phi - delta);

            // m = ba x cb
            const double m_x = ba_y*cb_z - ba_z*cb_y;
            const double m_y = ba_z*cb_x - ba_x*cb_z;
            const double m_z = ba_x*cb_y - ba_y*cb_x;

            // n = cb x dc

            const double n_x = cb_y*dc_z - cb_z*dc_y;
            const double n_y = cb_z*dc_x - cb_x*dc_z;
            const double n_z = cb_x*dc_y - cb_y*dc_x;

            const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
            if (mod_m_sq < 1e-24) continue;

            const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
            if (mod_n_sq < 1e-24) continue;

            const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
            if (mod_bc_sq < 1e-24) continue;

            const double mod_bc = std::sqrt(mod_bc_sq);

            const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
            const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

            const double fax = tau*mod_bc/mod_m_sq*m_x;
            const double fay = tau*mod_bc/mod_m_sq*m_y;
            const double faz = tau*mod_bc/mod_m_sq*m_z;

            const double fdx = -tau*mod_bc/mod_n_sq*n_x;
            const double fdy = -tau*mod_bc/mod_n_sq*n_y;
            const double fdz = -tau*mod_bc/mod_n_sq*n_z;

            const double fbx = (dot_ba_cb/mod_bc_sq - 1)*fax - dot_dc_cb/mod_bc_sq*fdx;
            const double fby = (dot_ba_cb/mod_bc_sq - 1)*fay - dot_dc_cb/mod_bc_sq*fdy;
            const double fbz = (dot_ba_cb/mod_bc_sq - 1)*faz - dot_dc_cb/mod_bc_sq*fdz;

            const double fcx = (dot_dc_cb/mod_bc_sq - 1)*fdx - dot_ba_cb/mod_bc_sq*fax;
            const double fcy = (dot_dc_cb/mod_bc_sq - 1)*fdy - dot_ba_cb/mod_bc_sq*fay;
            const double fcz = (dot_dc_cb/mod_bc_sq - 1)*fdz - dot_ba_cb/mod_bc_sq*faz;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
            forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
            forces_[3*atomDIndex] += fdx; forces_[3*atomDIndex+1] += fdy; forces_[3*atomDIndex+2] += fdz;
        }
    }

    void System::calculate_forces_harmonicImpropers() {
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (auto& dih: topology_.get_harmonic_impropers()) {
            const int atomAIndex = dih.get_atomA_index();
            const int atomBIndex = dih.get_atomB_index();
            const int atomCIndex = dih.get_atomC_index();
            const int atomDIndex = dih.get_atomD_index();

            const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                    3 * atomAIndex + 2];
            const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                    3 * atomBIndex + 2];
            const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                    3 * atomCIndex + 2];
            const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                    3 * atomDIndex + 2];

            // ba
            const double ba_x = x1 - x2;
            const double ba_y = y1 - y2;
            const double ba_z = z1 - z2;

            // cb
            const double cb_x = x2 - x3;
            const double cb_y = y2 - y3;
            const double cb_z = z2 - z3;

            // dc
            const double dc_x = x3 - x4;
            const double dc_y = y3 - y4;
            const double dc_z = z3 - z4;


            const double k = dih.get_IMP_force_constant();
            const double psi0 = dih.get_IMP_phase_value();
            const double psi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

            double delta = psi - psi0;
            while (delta <= -M_PI) delta += 2.0 * M_PI;
            while (delta > M_PI)   delta -= 2.0 * M_PI;
            const double tau = -2.0 * k * delta;

            // m = ba x cb
            const double m_x = ba_y*cb_z - ba_z*cb_y;
            const double m_y = ba_z*cb_x - ba_x*cb_z;
            const double m_z = ba_x*cb_y - ba_y*cb_x;

            // n = cb x dc

            const double n_x = cb_y*dc_z - cb_z*dc_y;
            const double n_y = cb_z*dc_x - cb_x*dc_z;
            const double n_z = cb_x*dc_y - cb_y*dc_x;

            const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
            if (mod_m_sq < 1e-24) continue;

            const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
            if (mod_n_sq < 1e-24) continue;

            const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
            if (mod_bc_sq < 1e-24) continue;

            const double mod_bc = std::sqrt(mod_bc_sq);

            const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
            const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

            const double fax = tau*mod_bc/mod_m_sq*m_x;
            const double fay = tau*mod_bc/mod_m_sq*m_y;
            const double faz = tau*mod_bc/mod_m_sq*m_z;

            const double fdx = -tau*mod_bc/mod_n_sq*n_x;
            const double fdy = -tau*mod_bc/mod_n_sq*n_y;
            const double fdz = -tau*mod_bc/mod_n_sq*n_z;

            const double fbx = (dot_ba_cb/mod_bc_sq - 1)*fax - dot_dc_cb/mod_bc_sq*fdx;
            const double fby = (dot_ba_cb/mod_bc_sq - 1)*fay - dot_dc_cb/mod_bc_sq*fdy;
            const double fbz = (dot_ba_cb/mod_bc_sq - 1)*faz - dot_dc_cb/mod_bc_sq*fdz;

            const double fcx = (dot_dc_cb/mod_bc_sq - 1)*fdx - dot_ba_cb/mod_bc_sq*fax;
            const double fcy = (dot_dc_cb/mod_bc_sq - 1)*fdy - dot_ba_cb/mod_bc_sq*fay;
            const double fcz = (dot_dc_cb/mod_bc_sq - 1)*fdz - dot_ba_cb/mod_bc_sq*faz;

            forces_[3*atomAIndex] += fax; forces_[3*atomAIndex+1] += fay; forces_[3*atomAIndex+2] += faz;
            forces_[3*atomBIndex] += fbx; forces_[3*atomBIndex+1] += fby; forces_[3*atomBIndex+2] += fbz;
            forces_[3*atomCIndex] += fcx; forces_[3*atomCIndex+1] += fcy; forces_[3*atomCIndex+2] += fcz;
            forces_[3*atomDIndex] += fdx; forces_[3*atomDIndex+1] += fdy; forces_[3*atomDIndex+2] += fdz;
        }
    }

    void System::calculate_forces_cmap()
    {
        const std::vector<double>& coordinates = topology_.get_coordinates();
        for (auto& cmap : topology_.get_cmaps()) {

            const int atomAIndex = cmap.get_atomA_index();
            const int atomBIndex = cmap.get_atomB_index();
            const int atomCIndex = cmap.get_atomC_index();
            const int atomDIndex = cmap.get_atomD_index();
            const int atomEIndex = cmap.get_atomE_index();

            const double x1 = coordinates[3 * atomAIndex], y1 = coordinates[3 * atomAIndex + 1], z1 = coordinates[
                    3 * atomAIndex + 2];
            const double x2 = coordinates[3 * atomBIndex], y2 = coordinates[3 * atomBIndex + 1], z2 = coordinates[
                    3 * atomBIndex + 2];
            const double x3 = coordinates[3 * atomCIndex], y3 = coordinates[3 * atomCIndex + 1], z3 = coordinates[
                    3 * atomCIndex + 2];
            const double x4 = coordinates[3 * atomDIndex], y4 = coordinates[3 * atomDIndex + 1], z4 = coordinates[
                    3 * atomDIndex + 2];
            const double x5 = coordinates[3 * atomEIndex], y5 = coordinates[3 * atomEIndex + 1], z5 = coordinates[
                    3 * atomEIndex + 2];


            const double phi = Metrics::dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);
            const double psi = Metrics::dihedral(x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5);

            // Calculate slopes (gradients)
            const int set_index = cmap.get_parameter_set();
            const std::vector<int>& resolutions = topology_.get_charmm_cmap_resolutions();
            const int resolution = resolutions[set_index - 1];
            const std::vector<double>& full_grid = topology_.get_cmap_grid_data(set_index);

            auto [fst, snd] = cmap.return_gradient_bicubic(phi, psi, resolution, topology_.get_Cmap_Coefficient_Matrix_bicubic_spline(), full_grid);
            const double dEdPhi = fst;
            const double dEdPsi = snd;



            const double ba_x = x1 - x2;
            const double ba_y = y1 - y2;
            const double ba_z = z1 - z2;

            // cb
            const double cb_x = x2 - x3;
            const double cb_y = y2 - y3;
            const double cb_z = z2 - z3;

            // dc
            const double dc_x = x3 - x4;
            const double dc_y = y3 - y4;
            const double dc_z = z3 - z4;

            //ed
            const double ed_x = x4 - x5;
            const double ed_y = y4 - y5;
            const double ed_z = z4 - z5;


            // m = ba x cb
            const double m_x = ba_y*cb_z - ba_z*cb_y;
            const double m_y = ba_z*cb_x - ba_x*cb_z;
            const double m_z = ba_x*cb_y - ba_y*cb_x;

            // n = cb x dc

            const double n_x = cb_y*dc_z - cb_z*dc_y;
            const double n_y = cb_z*dc_x - cb_x*dc_z;
            const double n_z = cb_x*dc_y - cb_y*dc_x;

            // o = dc x ed
            const double o_x = dc_y*ed_z - dc_z*ed_y;
            const double o_y = dc_z*ed_x - dc_x*ed_z;
            const double o_z = dc_x*ed_y - dc_y*ed_x;


            const double mod_m_sq = m_x*m_x + m_y*m_y + m_z*m_z;
            if (mod_m_sq < 1e-24) continue;

            const double mod_n_sq = n_x*n_x + n_y*n_y + n_z*n_z;
            if (mod_n_sq < 1e-24) continue;

            const double mod_o_sq = o_x*o_x + o_y*o_y + o_z*o_z;
            if (mod_o_sq < 1e-24) continue;

            const double mod_bc_sq = cb_x*cb_x + cb_y*cb_y + cb_z*cb_z;
            if (mod_bc_sq < 1e-24) continue;

            const double mod_cd_sq = dc_x*dc_x + dc_y*dc_y + dc_z*dc_z;
            if (mod_cd_sq < 1e-24) continue;

            const double mod_bc = std::sqrt(mod_bc_sq);
            const double mod_cd = std::sqrt(mod_cd_sq);

            const double dot_ba_cb = ba_x*cb_x + ba_y*cb_y + ba_z*cb_z;
            const double dot_dc_cb = dc_x*cb_x + dc_y*cb_y + dc_z*cb_z;

            const double dot_cb_dc = cb_x*dc_x + cb_y*dc_y + cb_z*dc_z;
            const double dot_ed_dc = ed_x*dc_x + ed_y*dc_y + ed_z*dc_z;




            const double fax_phi = -dEdPhi*mod_bc/mod_m_sq*m_x;
            const double fay_phi = -dEdPhi*mod_bc/mod_m_sq*m_y;
            const double faz_phi = -dEdPhi*mod_bc/mod_m_sq*m_z;

            const double fdx_phi = dEdPhi*mod_bc/mod_n_sq*n_x;
            const double fdy_phi = dEdPhi*mod_bc/mod_n_sq*n_y;
            const double fdz_phi = dEdPhi*mod_bc/mod_n_sq*n_z;

            const double fbx_phi = (dot_ba_cb/mod_bc_sq - 1)*fax_phi - dot_dc_cb/mod_bc_sq*fdx_phi;
            const double fby_phi = (dot_ba_cb/mod_bc_sq - 1)*fay_phi - dot_dc_cb/mod_bc_sq*fdy_phi;
            const double fbz_phi = (dot_ba_cb/mod_bc_sq - 1)*faz_phi - dot_dc_cb/mod_bc_sq*fdz_phi;

            const double fcx_phi = (dot_dc_cb/mod_bc_sq - 1)*fdx_phi - dot_ba_cb/mod_bc_sq*fax_phi;
            const double fcy_phi = (dot_dc_cb/mod_bc_sq - 1)*fdy_phi - dot_ba_cb/mod_bc_sq*fay_phi;
            const double fcz_phi = (dot_dc_cb/mod_bc_sq - 1)*fdz_phi - dot_ba_cb/mod_bc_sq*faz_phi;


            forces_[3*atomAIndex] += fax_phi; forces_[3*atomAIndex+1] += fay_phi; forces_[3*atomAIndex+2] += faz_phi;
            forces_[3*atomBIndex] += fbx_phi; forces_[3*atomBIndex+1] += fby_phi; forces_[3*atomBIndex+2] += fbz_phi;
            forces_[3*atomCIndex] += fcx_phi; forces_[3*atomCIndex+1] += fcy_phi; forces_[3*atomCIndex+2] += fcz_phi;
            forces_[3*atomDIndex] += fdx_phi; forces_[3*atomDIndex+1] += fdy_phi; forces_[3*atomDIndex+2] += fdz_phi;

            const double fbx_psi = -dEdPsi*mod_cd/mod_n_sq*n_x;
            const double fby_psi = -dEdPsi*mod_cd/mod_n_sq*n_y;
            const double fbz_psi = -dEdPsi*mod_cd/mod_n_sq*n_z;

            const double fex_psi = dEdPsi*mod_cd/mod_o_sq*o_x;
            const double fey_psi = dEdPsi*mod_cd/mod_o_sq*o_y;
            const double fez_psi = dEdPsi*mod_cd/mod_o_sq*o_z;

            const double fcx_psi = (dot_cb_dc/mod_cd_sq - 1)*fbx_psi - dot_ed_dc/mod_cd_sq*fex_psi;
            const double fcy_psi = (dot_cb_dc/mod_cd_sq - 1)*fby_psi - dot_ed_dc/mod_cd_sq*fey_psi;
            const double fcz_psi = (dot_cb_dc/mod_cd_sq - 1)*fbz_psi - dot_ed_dc/mod_cd_sq*fez_psi;

            const double fdx_psi = (dot_ed_dc/mod_cd_sq - 1)*fex_psi - dot_cb_dc/mod_cd_sq*fbx_psi;
            const double fdy_psi = (dot_ed_dc/mod_cd_sq - 1)*fey_psi - dot_cb_dc/mod_cd_sq*fby_psi;
            const double fdz_psi = (dot_ed_dc/mod_cd_sq - 1)*fez_psi - dot_cb_dc/mod_cd_sq*fbz_psi;

            forces_[3*atomBIndex] += fbx_psi; forces_[3*atomBIndex+1] += fby_psi; forces_[3*atomBIndex+2] += fbz_psi;
            forces_[3*atomCIndex] += fcx_psi; forces_[3*atomCIndex+1] += fcy_psi; forces_[3*atomCIndex+2] += fcz_psi;
            forces_[3*atomDIndex] += fdx_psi; forces_[3*atomDIndex+1] += fdy_psi; forces_[3*atomDIndex+2] += fdz_psi;
            forces_[3*atomEIndex] += fex_psi; forces_[3*atomEIndex+1] += fey_psi; forces_[3*atomEIndex+2] += fez_psi;
        }
    }
}