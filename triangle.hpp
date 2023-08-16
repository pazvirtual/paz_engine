#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include <array>

namespace paz
{
    class Triangle
    {
        const double x0, y0, z0;
        std::array<double, 3> basisX; // Negative normal vector.
        std::array<double, 3> basisY;
        std::array<double, 3> basisZ;
        double z1t, y2t, z2t;

    public:
        Triangle(double x0, double y0, double z0, double x1, double y1, double
            z1, double x2, double y2, double z2);
        double dist(double x, double y, double z) const;
        void collide(double x, double y, double z, double& nx, double& ny,
            double& nz, double& d) const;
    };
}

#endif
