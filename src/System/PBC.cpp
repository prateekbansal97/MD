//
// Created by Prateek Bansal on 12/20/25.
//

#include "System/System.h"

namespace md
{

    void System::set_ee_skin(const double eeskin)
    {
        ee_skin_ = eeskin;
        ee_list_cutoff_ = ee_skin_ + ee_cutoff_;
        ee_list_cutoff2_ = ee_list_cutoff_*ee_list_cutoff_;
    }

    int System::get_Cell_ID_lj(int ix, int iy, int iz) const {
        return ix + nCells_x_lj_*(iy + nCells_y_lj_*iz);
    }

    int System::get_Cell_ID_ee(int ix, int iy, int iz) const {
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

    double System::min_image_1d(double d, double L) {
        // move into [-L/2, L/2]
        d -= L * std::round(d / L);
        return d;
    }

    double  System::wrap(double x, double Lx) {
        x -= Lx * std::floor(x / Lx);
        return x;
    }


}