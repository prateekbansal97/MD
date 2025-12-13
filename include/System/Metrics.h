//
// Created by Prateek Bansal on 12/13/25.
//

#ifndef MD_METRICS_H
#define MD_METRICS_H

namespace Metrics
{
    double distance(double x1, double y1, double z1, double x2, double y2, double z2);
    double angle(double x1, double y1, double z1, double x2, double y2, double z2,  double x3, double y3, double z3);
    double dihedral(double x1, double y1, double z1, double x2, double y2, double z2,
                    double x3, double y3, double z3, double x4, double y4, double z4);
}
#endif //MD_METRICS_H