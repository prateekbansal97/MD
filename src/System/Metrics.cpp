//
// Created by Prateek Bansal on 12/13/25.
//

#include "../../include/System/Metrics.h"
#include <cmath>
#include <algorithm>

namespace Metrics
{

double distance(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
{
    const double dx = x1 - x2;
    const double dy = y1 - y2;
    const double dz = z1 - z2;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double angle(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,  const double x3, const double y3, const double z3)
{
    const double dx1 = x1 - x2;
    const double dy1 = y1 - y2;
    const double dz1 = z1 - z2;

    const double dx2 = x3 - x2;
    const double dy2 = y3 - y2;
    const double dz2 = z3 - z2;

    const double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
    const double mag1_sq = dx1*dx1 + dy1*dy1 + dz1*dz1;
    const double mag2_sq = dx2*dx2 + dy2*dy2 + dz2*dz2;
    if (constexpr double eps = 1e-14; mag1_sq < eps || mag2_sq < eps) {
        return std::numeric_limits<double>::quiet_NaN(); // or 0.0, or throw
    }
    const double inv_mag = 1.0 / std::sqrt(mag1_sq * mag2_sq);
    double c = dot * inv_mag;
    c = std::clamp(c, -1.0, 1.0);
    return std::acos(c);
}

double dihedral(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2,
                const double x3, const double y3, const double z3, const double x4, const double y4, const double z4)
{
    // Calculating P1, P2, P3 as 2 - 1, 3 - 2, 4 - 3
    const double dx1 = x2 - x1;
    const double dy1 = y2 - y1;
    const double dz1 = z2 - z1;
    const double dx2 = x3 - x2;
    const double dy2 = y3 - y2;
    const double dz2 = z3 - z2;

    const double dx3 = x4 - x3;
    const double dy3 = y4 - y3;
    const double dz3 = z4 - z3;

    // Calculating n1 = P1 x P2
    const double n1x = dy1*dz2 - dz1*dy2;
    const double n1y = dx2*dz1 - dx1*dz2; //-1 * (dx1*dz2 - dx2*dz1);
    const double n1z = dx1*dy2 - dy1*dx2;

    // Calculating n2 = P2 x P3
    const double n2x = dy2*dz3 - dz2*dy3;
    const double n2y = dx3*dz2 - dx2*dz3;//-1 * (dx2*dz3 - dx3*dz2);
    const double n2z = dx2*dy3 - dy2*dx3;

    const double magn1 = n1x*n1x + n1y*n1y + n1z*n1z;
    const double magn2 = n2x*n2x + n2y*n2y + n2z*n2z;
    const double magb1 = std::sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

    constexpr double eps2 = 1e-28;

    if (constexpr double eps = 1e-14; magn1 < eps2 || magn2 < eps2 || magb1 < eps) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double x = n1x*n2x + n1y*n2y + n1z*n2z; // n1 dot n2

    const double val = dx1*n2x + dy1*n2y + dz1*n2z; // n2 dot b1
    const double y = val*magb1;
    return std::atan2(y, x);
}
}
