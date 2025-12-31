//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"
#include "/usr/local/Cellar/fftw/3.3.10_2/include/fftw3.h"
#include <array>
#include <complex>

namespace md
{
    void System::set_ee_skin(const double eeskin)
    {
        ee_skin_ = eeskin;
        ee_list_cutoff_ = ee_skin_ + ee_cutoff_;
        ee_list_cutoff2_ = ee_list_cutoff_*ee_list_cutoff_;
    }

    int System::get_Cell_ID_lj(const int ix, const int iy, const int iz) const {
        return ix + nCells_x_lj_*(iy + nCells_y_lj_*iz);
    }

    int System::get_Cell_ID_ee(const int ix, const int iy, const int iz) const {
        return ix + nCells_x_ee_*(iy + nCells_y_ee_*iz);
    }

    void System::set_box(const double Lx, const double Ly, const double Lz) {
        boxLx_ = Lx; boxLy_ = Ly; boxLz_ = Lz;
        pbc_enabled_ = (Lx > 0 && Ly > 0 && Lz > 0);
    }

    int System::cell_number_along_dim(const double x, const double L, const int nCelldim)
    {
        // x is assumed wrapped into [0, L)
        int ix = static_cast<int>(x * nCelldim / L);   // 0..nCells-1
        if (ix >= nCelldim) ix = nCelldim - 1;           // safety for x ~ L due to fp
        if (ix < 0) ix = 0;
        return ix;
    }

    void System::init_box()
    {
        if (!pbc_enabled_) return;
        nCells_x_lj_ = std::max(1, static_cast<int>(std::floor(boxLx_ / get_lj_list_cutoff())));
        nCells_y_lj_ = std::max(1, static_cast<int>(std::floor(boxLy_ / get_lj_list_cutoff())));
        nCells_z_lj_ = std::max(1, static_cast<int>(std::floor(boxLz_ / get_lj_list_cutoff())));

        lcell_x_lj_ = boxLx_ / nCells_x_lj_;
        lcell_y_lj_ = boxLy_ / nCells_y_lj_;
        lcell_z_lj_ = boxLz_ / nCells_z_lj_;

        nCells_lj_ = nCells_x_lj_ * nCells_y_lj_ * nCells_z_lj_;

        nCells_x_ee_ = std::max(1, static_cast<int>(std::floor(boxLx_ / get_ee_list_cutoff())));
        nCells_y_ee_ = std::max(1, static_cast<int>(std::floor(boxLy_ / get_ee_list_cutoff())));
        nCells_z_ee_ = std::max(1, static_cast<int>(std::floor(boxLz_ / get_ee_list_cutoff())));

        lcell_x_ee_ = boxLx_ / nCells_x_ee_;
        lcell_y_ee_ = boxLy_ / nCells_y_ee_;
        lcell_z_ee_ = boxLz_ / nCells_z_ee_;

        nCells_ee_ = nCells_x_ee_ * nCells_y_ee_ * nCells_z_ee_;
    }

    void System::apply_min_image(double& dx, double& dy, double& dz) const {
        if (!pbc_enabled_) return;
        dx = min_image_1d(dx, boxLx_);
        dy = min_image_1d(dy, boxLy_);
        dz = min_image_1d(dz, boxLz_);
    }

    double System::min_image_1d(double d, const double L) {
        // move into [-L/2, L/2]
        d -= L * std::round(d / L);
        return d;
    }

    double  System::wrap(double x, const double Lx) {
        x -= Lx * std::floor(x / Lx);
        return x;
    }

    void System::set_ewald_error_tolerance(const double tolerance)
    {
        ewald_error_tolerance_ = tolerance;
    }

    void System::calculate_ewald_alpha()
    {
        ewald_alpha_ = std::sqrt(std::log( 2 * ewald_error_tolerance_ )*-1)/ee_cutoff_;
    }

    void System::calculate_ewald_nnodes()
    {
        newald_mesh_x = static_cast<int>((2*ewald_alpha_*boxLx_)/(3*std::pow(ewald_error_tolerance_, 0.2)));
        newald_mesh_y = static_cast<int>((2*ewald_alpha_*boxLy_)/(3*std::pow(ewald_error_tolerance_, 0.2)));
        newald_mesh_z = static_cast<int>((2*ewald_alpha_*boxLz_)/(3*std::pow(ewald_error_tolerance_, 0.2)));
        newald_total = newald_mesh_x*newald_mesh_y*newald_mesh_z;
    }

    double System::calculate_center_of_mass_along_direction(const std::string_view direction)
    {
        if (direction.empty()) {
            std::cerr << "ERROR: empty direction\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        const std::vector<double>& coordinates = topology_.get_coordinates();
        const char axis = static_cast<char>(std::toupper(static_cast<unsigned char>(direction[0])));
        int offset = -1;
        if (axis == 'X') { offset = 0; }
        else if (axis == 'Y') { offset = 1; }
        else if (axis == 'Z') { offset = 2; }
        else {std::cerr << "ERROR: unknown direction \"" << direction << "\"" << std::endl;}

        double total_mass = 0;
        double weighted_sum = 0;

        for (size_t i = 0; i < natoms_; ++i)
        {
            const double mass = topology_.get_atom_list()[i].get_mass();
            const double coordinate = coordinates[i*3 + offset];

            weighted_sum +=  coordinate * mass;
            total_mass += mass;
        }

        if (total_mass == 0.0) {
            std::cerr << "ERROR: total mass is zero\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        return weighted_sum/total_mass;
    }

    void System::move_box(const std::string_view direction, const double displacement)
    {
        std::vector<double> coordinates = topology_.get_coordinates();
        const char axis = static_cast<char>(std::toupper(static_cast<unsigned char>(direction[0])));
        int offset = -1;
        if (axis == 'X') { offset = 0; }
        else if (axis == 'Y') { offset = 1; }
        else if (axis == 'Z') { offset = 2; }

        for (size_t i = 0; i < natoms_; ++i)
        {
            coordinates[i*3 + offset] -= displacement;
        }
    }

    void System::calculate_pme_grid_point_positions()
    {
        const double low_x = -1 * boxLx_ / 2;
        const double low_y = -1 * boxLy_ / 2;
        const double low_z = -1 * boxLz_ / 2;

        for (int i = 0; i < newald_mesh_x; ++i)
        {
            double node_position_x = low_x + i * boxLx_/newald_mesh_x;
            pme_grid_point_positions_x_.push_back(node_position_x);
        }

        for (int i = 0; i < newald_mesh_y; ++i)
        {
            double node_position_y = low_y + i * boxLy_/newald_mesh_y;
            pme_grid_point_positions_y_.push_back(node_position_y);
        }

        for (int i = 0; i < newald_mesh_z; ++i)
        {
            double node_position_z = low_z + i * boxLz_/newald_mesh_z;
            pme_grid_point_positions_z_.push_back(node_position_z);
        }
        // pme_grid_charge_x_.assign(newald_mesh_x, 0.0);
        // pme_grid_charge_y_.assign(newald_mesh_y, 0.0);
        // pme_grid_charge_z_.assign(newald_mesh_z, 0.0);
        // for (const auto i: pme_grid_point_positions_x_)
        // {
        //     std::cout << i << " ";
        // }
    }

    void System::spread_charge_pme()
    {

        const double scale_x = newald_mesh_x / boxLx_;
        const double scale_y = newald_mesh_y / boxLy_;
        const double scale_z = newald_mesh_z / boxLz_;

        pme_grid_charge_.assign(newald_total, 0.0);
        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (size_t i = 0; i < natoms_; ++i)
        {
            const double charge = topology_.get_atom_list()[i].get_partial_charge();
            if (charge == 0.0) continue;

            const double x = wrap(coordinates[3*i], boxLx_);
            const double y = wrap(coordinates[3*i + 1], boxLy_);
            const double z = wrap(coordinates[3*i + 2], boxLz_);

            const double u_x = x * scale_x;
            const double u_y = y * scale_y;
            const double u_z = z * scale_z;

            const int idx_x = static_cast<int>(std::floor(u_x));
            const int idx_y = static_cast<int>(std::floor(u_y));
            const int idx_z = static_cast<int>(std::floor(u_z));

            const double dx = u_x - idx_x;
            const double dy = u_y - idx_y;
            const double dz = u_z - idx_z;

            auto theta_x = calculate_b_spline_coeffs_order5(dx);
            auto theta_y = calculate_b_spline_coeffs_order5(dy);
            auto theta_z = calculate_b_spline_coeffs_order5(dz);

            constexpr int order = 5;

            for (int ix = 0; ix < order; ++ix)
            {
                constexpr int offset = -2;

                int kx = idx_x + offset + ix;
                kx = (kx % newald_mesh_x + newald_mesh_x) % newald_mesh_x;

                const double wx = theta_x[ix]; // Optimization: Load x-weight once

                for (int iy = 0; iy < order; ++iy)
                {
                    int ky = idx_y + offset + iy;
                    ky = (ky % newald_mesh_y + newald_mesh_y) % newald_mesh_y;

                    const double wxy = wx * theta_y[iy];

                    for (int iz = 0; iz < order; ++iz)
                    {
                        int kz = idx_z + offset + iz;
                        kz = (kz % newald_mesh_z + newald_mesh_z) % newald_mesh_z;

                        const int index = kx * (newald_mesh_y * newald_mesh_z) +
                                    ky * (newald_mesh_z) +
                                    kz;

                        pme_grid_charge_[index] += charge * wxy * theta_z[iz];
                    }
                }
            }
        }
    }

    double System::findClosestSorted(const std::vector<double>& vec, const double value) {
        if (vec.empty()) {
            throw std::runtime_error("Vector cannot be empty");
        }

        const auto it_upper = std::ranges::lower_bound(vec, value);

        if (it_upper == vec.begin()) {
            return *it_upper; // Value is smaller than the smallest element
        }
        if (it_upper == vec.end()) {
            return *std::prev(it_upper); // Value is larger than the largest element
        }

        if (const auto it_lower = std::prev(it_upper); std::abs(*it_upper - value) < std::abs(*it_lower - value)) {
            return *it_upper;
        } else {
            return *it_lower;
        }
    }

    size_t System::findClosestIndexSorted(const std::vector<double>& vec, const double value) {
        if (vec.empty()) {
            throw std::runtime_error("Vector cannot be empty");
        }

        const auto it = std::ranges::lower_bound(vec, value);

        // Case 1: value is smaller than or equal to all elements
        if (it == vec.begin()) {
            return 0;
        }

        // Case 2: value is larger than all elements
        if (it == vec.end()) {
            return vec.size() - 1;
        }

        // Case 3: value is between (it - 1) and it
        const double diff1 = std::abs(value - *it);

        if (const double diff2 = std::abs(value - *(it - 1)); diff2 < diff1) {
            return static_cast<size_t>(std::distance(vec.begin(), it - 1));
        }

        return static_cast<size_t>(std::distance(vec.begin(), it));
    }


    // std::vector<double> System::calculate_b_spline_coeffs(const double diff)
    // {
    //     std::vector coeffs_ = {diff, 1-diff};
    //
    //     for (int k = 3; k <= 5; k++)
    //     {
    //         std::vector<double> new_coeffs_(k, 0.0);
    //         const double divterm = 1.0 / (k - 1);
    //
    //         for (int i = 0; i < k; i++)
    //         {
    //             const double H_i = (i < coeffs_.size()) ? coeffs_[i] : 0.0;
    //             const double H_im1 = (i - 1 >= 0) ? coeffs_[i - 1] : 0.0;
    //
    //             new_coeffs_[i] = divterm * (diff * H_i + (k - diff)*H_im1);
    //         }
    //         coeffs_ = new_coeffs_;
    //     }
    //
    //     return coeffs_;
    // }

    // #include <array>

    // Calculates B-Spline coefficients for PME (Order 5)
    // Input 'u': The fractional part of the coordinate (x - floor(x)), range [0, 1)
    // Returns: 5 weights corresponding to grid points [i-2, i-1, i, i+1, i+2] relative to floor

    std::array<double, 5> System::calculate_b_spline_coeffs_order5(const double u)
    {
        std::array<double, 5> theta{};

        // Current Order 2 weights (Linear)
        const double c0 = 1.0 - u;
        const double c1 = u;

        // Iteration k=3
        const double w0 = (1.0/2.0) * ((1.0-u)*c0);
        const double w1 = (1.0/2.0) * ((1.0+u)*c0 + (2.0-u)*c1);
        const double w2 = (1.0/2.0) * (             (1.0+u)*c1);

        // Iteration k=4
        const double x0 = (1.0/3.0) * ((1.0-u)*w0);
        const double x1 = (1.0/3.0) * ((2.0+u)*w0 + (2.0-u)*w1);
        const double x2 = (1.0/3.0) * ((1.0+u)*w1 + (3.0-u)*w2);
        const double x3 = (1.0/3.0) * (             (2.0+u)*w2);

        // Iteration k=5 (Final)
        theta[0] = (1.0/4.0) * ((1.0-u)*x0);
        theta[1] = (1.0/4.0) * ((3.0+u)*x0 + (2.0-u)*x1);
        theta[2] = (1.0/4.0) * ((2.0+u)*x1 + (3.0-u)*x2);
        theta[3] = (1.0/4.0) * ((1.0+u)*x2 + (4.0-u)*x3);
        theta[4] = (1.0/4.0) * (             (3.0+u)*x3);

        return theta;
    }


    std::pair<std::array<double, 5>, std::array<double, 5>> System::calculate_b_spline_coeffs_and_derivs_order5(const double u)
    {
        std::array<double, 5> theta{};
        std::array<double, 5> dtheta{};

        // --- Order 2 (Linear) ---
        const double w2_0 = 1.0 - u;
        const double w2_1 = u;

        // --- Order 3 (Quadratic) ---
        const double w3_0 = 0.5 * (1.0 - u) * w2_0;
        const double w3_1 = 0.5 * ((1.0 + u) * w2_0 + (2.0 - u) * w2_1);
        const double w3_2 = 0.5 * (1.0 + u) * w2_1;

        // --- Order 4 (Cubic) ---
        constexpr double div4 = 1.0 / 3.0;
        const double w4_0 = div4 * (1.0 - u) * w3_0;
        const double w4_1 = div4 * ((2.0 + u) * w3_0 + (2.0 - u) * w3_1);
        const double w4_2 = div4 * ((1.0 + u) * w3_1 + (3.0 - u) * w3_2);
        const double w4_3 = div4 * (2.0 + u) * w3_2;

        // --- Order 5 (Quintic) Values ---
        constexpr double div5 = 0.25;
        theta[0] = div5 * (1.0 - u) * w4_0;
        theta[1] = div5 * ((3.0 + u) * w4_0 + (2.0 - u) * w4_1);
        theta[2] = div5 * ((2.0 + u) * w4_1 + (3.0 - u) * w4_2);
        theta[3] = div5 * ((1.0 + u) * w4_2 + (4.0 - u) * w4_3);
        theta[4] = div5 * (3.0 + u) * w4_3;

        // --- Order 5 Derivatives ---
        // dM_5(u)/du = M_4(u) - M_4(u-1)
        // Note: The w4 variables calculated above are M_4 shifted.
        // Standard analytic derivatives for cardinality (Order 5) centered at 0:
        // This is equivalent to M4(u + shift) logic.
        // The simplest robust way without re-deriving the recurrence is:
        dtheta[0] = -w4_0;
        dtheta[1] = w4_0 - w4_1;
        dtheta[2] = w4_1 - w4_2;
        dtheta[3] = w4_2 - w4_3;
        dtheta[4] = w4_3;

        return {theta, dtheta};
    }


    void System::init_pme_fft()
    {
        if (fft_initialized_) return;
        const size_t complex_size = newald_mesh_x * newald_mesh_y * (newald_mesh_z / 2 + 1);

        pme_grid_complex_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * complex_size));

        if (!pme_grid_complex_) {
            throw std::runtime_error("Failed to allocate memory for FFTW complex grid.");
        }


        fft_forward_plan_ = fftw_plan_dft_r2c_3d(
        newald_mesh_x, newald_mesh_y, newald_mesh_z,
        pme_grid_charge_.data(),    // Input (Real)
        pme_grid_complex_,          // Output (Complex)
        FFTW_ESTIMATE
    );

        fft_backward_plan_ = fftw_plan_dft_c2r_3d(
        newald_mesh_x, newald_mesh_y, newald_mesh_z,
        pme_grid_complex_,          // Input (Complex)
        pme_grid_charge_.data(),    // Output (Real) - recycles the charge vector
        FFTW_ESTIMATE
    );

        fft_initialized_ = true;
        std::cout << "FFTW3 Plans Initialized." << std::endl;

    }

    void System::perform_forward_fft()
    {
        if (!fft_initialized_) {
            init_pme_fft();
        }

        // Execute the Real-to-Complex transform
        // Reads from pme_grid_charge_, writes to pme_grid_complex_
        fftw_execute(fft_forward_plan_);
    }

    System::~System()
    {
        if (fft_initialized_) {
            fftw_destroy_plan(fft_forward_plan_);
            fftw_destroy_plan(fft_backward_plan_);
            fftw_free(pme_grid_complex_);
        }
    }

    std::vector<double> System::compute_moduli(const int mesh_size)
    {

        std::vector<double> real_data(mesh_size, 0.0);
        const auto theta = calculate_b_spline_coeffs_order5(0.0);

        real_data[0] = theta[2]; // Center (offset 0)
        real_data[1] = theta[3]; // Offset +1
        real_data[2] = theta[4]; // Offset +2
        real_data[mesh_size - 1] = theta[1]; // Offset -1
        real_data[mesh_size - 2] = theta[0]; // Offset -2

        const int complex_size = mesh_size / 2 + 1;
        auto* complex_out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * complex_size));

        fftw_plan plan = fftw_plan_dft_r2c_1d(mesh_size, real_data.data(), complex_out, FFTW_ESTIMATE);
        fftw_execute(plan);

        std::vector<double> moduli(mesh_size); // Store for full range 0..N-1 for easier lookup later

        for (int i = 0; i < complex_size; ++i) {
            // Compute magnitude squared: Re^2 + Im^2
            const double mag2 = complex_out[i][0] * complex_out[i][0] +
                          complex_out[i][1] * complex_out[i][1];

            moduli[i] = mag2;

            if (i > 0 && i < (mesh_size + 1) / 2) {
                moduli[mesh_size - i] = mag2;
            }
        }

        fftw_destroy_plan(plan);
        fftw_free(complex_out);
        return moduli;
    }

    void System::init_pme_resources()
    {
        init_pme_fft();
        // Precompute B-spline moduli for all dimensions
        bspline_moduli_x_ = compute_moduli(newald_mesh_x);
        bspline_moduli_y_ = compute_moduli(newald_mesh_y);
        bspline_moduli_z_ = compute_moduli(newald_mesh_z);
    }

    void System::solve_poisson_kspace()
    {
        const double volume = boxLx_ * boxLy_ * boxLz_;
        const double factor_pre = 4.0 * M_PI / volume; // 4*pi*V
        const double beta_sq_inv = -1.0 / (4.0 * ewald_alpha_ * ewald_alpha_);
        const double fft_norm = 1.0 / static_cast<double>(newald_total);

        // Iterate over the complex grid (Nx, Ny, Nz/2 + 1)
        // Note: FFTW R2C layout stores only the non-redundant half of the last dimension.

        for (int i = 0; i < newald_mesh_x; ++i)
        {
            // Frequency kx: if i <= N/2, k=i; else k = i - N
            const int mx = (i <= newald_mesh_x / 2) ? i : i - newald_mesh_x;
            const double kx = 2.0 * M_PI * mx / boxLx_;

            for (int j = 0; j < newald_mesh_y; ++j)
            {
                const int my = (j <= newald_mesh_y / 2) ? j : j - newald_mesh_y;
                const double ky = 2.0 * M_PI * my / boxLy_;

                for (int k = 0; k < newald_mesh_z / 2 + 1; ++k)
                {
                    // Last dim is always 0..N/2
                    const double kz = 2.0 * M_PI * k / boxLz_;

                    // Calculate k^2
                    const double k_sq = kx*kx + ky*ky + kz*kz;

                    // 1. Handle k=0 case (Background charge neutralization)
                    if (k_sq < 1e-10) {
                        const int idx = i * (newald_mesh_y * (newald_mesh_z / 2 + 1)) +
                                  j * (newald_mesh_z / 2 + 1) +
                                  k;
                        pme_grid_complex_[idx][0] = 0.0;
                        pme_grid_complex_[idx][1] = 0.0;
                        continue;
                    }

                    // 2. Green's Function: exp(-k^2 / 4a^2) / k^2
                    const double influence_function = factor_pre * std::exp(k_sq * beta_sq_inv) / k_sq;

                    // 3. B-spline Moduli Correction (Inverse)
                    // We divide by |M|^2 for each dimension.
                    // Note: Since we use B-splines for both spread and force interpolation,
                    // the total correction is 1 / |M|^2. (Sometimes written as 1/|b|^2)

                    // Be careful with indices into moduli arrays!
                    // moduli arrays are size N, so we can use i, j directly (they wrap naturally).
                    // For k (z-dim), we are in 0..N/2, which maps 1:1.
                    const double b_corr = 1.0 / (bspline_moduli_x_[i] * bspline_moduli_y_[j] * bspline_moduli_z_[k]);

                    const double total_scale = influence_function * b_corr * fft_norm;
                    // 4. Apply to Grid
                    const int idx = i * (newald_mesh_y * (newald_mesh_z / 2 + 1)) +
                              j * (newald_mesh_z / 2 + 1) +
                              k;

                    pme_grid_complex_[idx][0] *= total_scale; // Real part
                    pme_grid_complex_[idx][1] *= total_scale; // Imag part
                }
            }
        }
    }

    void System::perform_backward_fft()
    {
        if (!fft_initialized_) {
            throw std::runtime_error("FFT not initialized before backward transform.");
        }
        // Complex -> Real
        // Result is stored in pme_grid_charge_ (which now holds Potential, not Charge)
        fftw_execute(fft_backward_plan_);
    }

    void System::gather_forces_pme()
    {
        const double scale_x = newald_mesh_x / boxLx_;
        const double scale_y = newald_mesh_y / boxLy_;
        const double scale_z = newald_mesh_z / boxLz_;

        const std::vector<double>& coordinates = topology_.get_coordinates();

        for (size_t i = 0; i < natoms_; ++i)
        {
            const double charge = topology_.get_atom_list()[i].get_partial_charge();
            if (charge == 0.0) continue;

            const double x = wrap(coordinates[3*i], boxLx_);
            const double y = wrap(coordinates[3*i + 1], boxLy_);
            const double z = wrap(coordinates[3*i + 2], boxLz_);

            const double u_x = x * scale_x;
            const double u_y = y * scale_y;
            const double u_z = z * scale_z;

            const int idx_x = static_cast<int>(std::floor(u_x));
            const int idx_y = static_cast<int>(std::floor(u_y));
            const int idx_z = static_cast<int>(std::floor(u_z));

            const double dx = u_x - idx_x;
            const double dy = u_y - idx_y;
            const double dz = u_z - idx_z;

            // Get Values AND Derivatives
            auto [theta_x, dtheta_x] = calculate_b_spline_coeffs_and_derivs_order5(dx);
            auto [theta_y, dtheta_y] = calculate_b_spline_coeffs_and_derivs_order5(dy);
            auto [theta_z, dtheta_z] = calculate_b_spline_coeffs_and_derivs_order5(dz);

            constexpr int order = 5;
            constexpr int offset = -2;

            double fx = 0.0;
            double fy = 0.0;
            double fz = 0.0;

            for (int ix = 0; ix < order; ++ix)
            {
                int kx = (idx_x + offset + ix) % newald_mesh_x;
                if (kx < 0) kx += newald_mesh_x;

                for (int iy = 0; iy < order; ++iy)
                {
                    int ky = (idx_y + offset + iy) % newald_mesh_y;
                    if (ky < 0) ky += newald_mesh_y;

                    for (int iz = 0; iz < order; ++iz)
                    {
                        int kz = (idx_z + offset + iz) % newald_mesh_z;
                        if (kz < 0) kz += newald_mesh_z;

                        const int index = kx * (newald_mesh_y * newald_mesh_z) +
                                          ky * (newald_mesh_z) +
                                          kz;

                        // The potential at this grid point
                        const double potential = pme_grid_charge_[index];

                        // F = -q * grad(Phi)
                        // dPhi/dx = potential * dMx/dx * My * Mz * scale_x

                        fx -= charge * potential * dtheta_x[ix] * theta_y[iy] * theta_z[iz] * scale_x;
                        fy -= charge * potential * theta_x[ix] * dtheta_y[iy] * theta_z[iz] * scale_y;
                        fz -= charge * potential * theta_x[ix] * theta_y[iy] * dtheta_z[iz] * scale_z;
                    }
                }
            }

            forces_[3*i + 0] += fx;
            forces_[3*i + 1] += fy;
            forces_[3*i + 2] += fz;
        }
    }
}

