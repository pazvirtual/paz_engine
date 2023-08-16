#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include <array>

namespace paz
{
    class Triangle
    {
        const double x0, y0, z0;
        const bool _degenerate;
        std::array<double, 3> basisX; // Negative normal vector.
        std::array<double, 3> basisY;
        std::array<double, 3> basisZ;
        double z1t, y2t, z2t;
        const std::array<double, 3> _centroid;
        double _radius;
        double dist_transformed(double xt, double yt, double zt, double&
            deltaYt, double& deltaZt) const;

    public:
        Triangle(double x0, double y0, double z0, double x1, double y1, double
            z1, double x2, double y2, double z2);
        void collide(double x, double y, double z, double radius, double& nx,
            double& ny, double& nz, double& d) const;
    };
}

#endif
