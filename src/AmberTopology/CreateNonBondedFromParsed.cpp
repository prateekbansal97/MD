//
// Created by Prateek Bansal on 12/21/25.
//

#include "AmberTopology/AmberTopology.h"

namespace md::AmberTopology
{
    void AmberTopology::create_excluded_atoms_list()
    {
        size_t start = 0;

        // Safety check: Ensure the list is large enough
        // This prevents a crash if the file is truncated
        const size_t total_expected = excluded_atoms_list_.size();

        for (auto& atom: atom_list_)
        {
            const unsigned long int excluded = atom.get_nExcluded_Atoms();

            // Bounds check
            if (start + excluded > total_expected) {
                std::cerr << "Error: Excluded atoms list is shorter than expected!" << std::endl;
                break;
            }

            for (size_t p = start; p < start + excluded; p++)
            {
                const int amber_index = excluded_atoms_list_[p];

                // AMBER index 0 is a placeholder/null.
                // If we find 0, we shouldn't add -1 to the list.
                auto excluded_vector = atom.get_excluded_atoms();
                if (amber_index > 0) {
                    excluded_vector.push_back(amber_index - 1);
                }
            }
            start += excluded;

            std::ranges::sort(atom.get_excluded_atoms());
            atom.get_excluded_atoms().erase(std::ranges::unique(atom.get_excluded_atoms()).begin(),
                                          atom.get_excluded_atoms().end());
        }
    }

    void AmberTopology::create_Cmap_Coefficient_Matrix_bicubic_spline()
    {
        static const double coeff_matrix[16][16] =
        {{0,   0,   0,   0,   0,  16,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0},
         {0,   0,   0,   0,  -8,   0,   8,   0,   0,   0,   0,   0,   0,   0,   0,   0},
         {0,   0,   0,   0,  16, -40,  32,  -8,   0,   0,   0,   0,   0,   0,   0,   0},
         {0,   0,   0,   0,  -8,  24, -24,   8,   0,   0,   0,   0,   0,   0,   0,   0},
         {0,  -8,   0,   0,   0,   0,   0,   0,   0,   8,   0,   0,   0,   0,   0,   0},
         {0,  -4,   0,   0,  -4,   4,   0,   0,   0,   0,   4,   0,   0,   0,   0,   0},
         {0,  32, -20,   0,   8,  -4,  -4,   0,   0, -24,  16,  -4,   0,   0,   0,   0},
         {0, -20,  12,   0,  -4,   0,   4,   0,   0,  16, -12,   4,   0,   0,   0,   0},
         {0,  16,   0,   0,   0, -40,   0,   0,   0,  32,   0,   0,   0,  -8,   0,   0},
         {0,   8,   0,   0,  32,  -4, -24,   0, -20,  -4,   0,   0,   0,   0,  -4,   0},
         {0, -64,  40,   0, -64,  96, -68,  24,  40, -68, -16, -16,   0,  24, -16,   0},
         {0,  40, -24,   0,  32, -52,  52, -24, -20,  40,  16,  16,   0, -16,  12,   0},
         {0,  -8,   0,   0,   0,  24,   0,   0,   0, -24,   0,   0,   0,   8,   0,   4},
         {0,  -4,   0,   0, -20,   0,  16,   0,  12,   4, -12,   0,   0,   0,   4,  -4},
         {0,  32, -20,   0,  40, -52,  40, -16, -24,  52, -52,  12,   0, -24,  16,  -4},
         {0, -20,  12,   0, -20,  28, -32,  16,  12, -32,  40, -12,   0,  16, -12,   4}};

        coeff_matrix_24x24_.clear();
        coeff_matrix_24x24_.reserve(5*24*24*16);

        for (int cmap_no = 0; cmap_no < 5; cmap_no++)
        {
            const int resolution = get_charmm_cmap_resolutions()[cmap_no];
            std::vector<double> grid = get_cmap_grid_data(cmap_no + 1);


            for (int row_grid = 0; row_grid < resolution; row_grid++)
            {
                for (int col_grid = 0; col_grid < resolution; col_grid++)
                {
                    int sixteen_indices[16][2] = {
                        {row_grid-1, col_grid-1}, {row_grid, col_grid-1}, {row_grid+1, col_grid-1}, {row_grid+2, col_grid-1},
                        {row_grid-1, col_grid},   {row_grid, col_grid},   {row_grid+1, col_grid},   {row_grid+2, col_grid},
                        {row_grid-1, col_grid+1}, {row_grid, col_grid+1}, {row_grid+1, col_grid+1}, {row_grid+2, col_grid+1},
                        {row_grid-1, col_grid+2}, {row_grid, col_grid+2}, {row_grid+1, col_grid+2}, {row_grid+2, col_grid+2}
                    };

                    for (auto & sixteen_index : sixteen_indices)
                    {
                        if (const int index_x = sixteen_index[0]; index_x < 0)  sixteen_index[0] += resolution;
                        else if (index_x >= resolution) sixteen_index[0] -= resolution;

                        if (const int index_y = sixteen_index[1]; index_y < 0) sixteen_index[1] += resolution;
                        else if (index_y >= resolution) sixteen_index[1] -= resolution;
                    }


                    std::vector<double> sixteen_energies(16, 0.0);
                    int ct = 0;
                    for (const auto& index: sixteen_indices)
                    {
                        const int index_x = index[0];
                        const int index_y = index[1];
                        sixteen_energies[ct] = grid[index_x*resolution + index_y];
                        ct++;
                    }

                    std::vector<double> coeffs(16, 0.0);
                    for (int row = 0; row < 16; ++row) {
                        // Iterate through each column (element) in the row
                        for (int col = 0; col < 16; ++col) {
                            coeffs[row] += coeff_matrix[row][col] * sixteen_energies[col] / 16;
                        }
                    }
                    coeff_matrix_24x24_.insert(coeff_matrix_24x24_.end(), coeffs.begin(), coeffs.end());
                }
            }
        }
    }

    static inline uint64_t pair_key(int a, int b) {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) | static_cast<uint32_t>(b);
    }

    void AmberTopology::build_lj14_pairlist() {
        lj14_pairs_.clear();
        lj14_pairs_.reserve(CosineDihedral_list_.size());

        scee14_scale_map_.clear();
        scee14_scale_map_.reserve(CosineDihedral_list_.size());

        for (const auto& dih : CosineDihedral_list_) {
            if (dih.get_improper()) continue;      // optional: usually skip impropers for 1-4
            if (dih.get_exclude_14()) continue;       // AMBER says "no 1-4 for this one"

            const int a = dih.get_atomA_index();
            const int d = dih.get_atomD_index();
            lj14_pairs_.insert(pair_key(a, d));

            const int t = dih.get_type();   // <-- rename to whatever you actually have

            double scale = 1.0;

            if (!scee_scale_factors_.empty() && t >= 0 && static_cast<size_t>(t) < scee_scale_factors_.size()) {
                scale = scee_scale_factors_[static_cast<size_t>(t)];
            }

            const uint64_t key = pair_key(a, d);
            if (auto it = scee14_scale_map_.find(key); it == scee14_scale_map_.end()) {
                scee14_scale_map_.emplace(key, scale);
            } else
            {
                // If multiple dihedrals generate the same 1-4 pair:
                // keep it only if consistent; otherwise prefer the *smaller* scaling (more conservative)
                if (const double prev = it->second; std::abs(prev - scale) > 1e-12) {
                    it->second = std::max(prev, scale);
                }
            }
        }
    }

    bool AmberTopology::is_14_pair(const int i, const int j) const {
        return lj14_pairs_.contains(pair_key(i, j));
    }

    double AmberTopology::get_scee_scale_for_pair(const int i, const int j) const {
        const auto it = scee14_scale_map_.find(pair_key(i, j));
        if (it == scee14_scale_map_.end()) return 1.0;   // default: unscaled
        return it->second;
    }
}