#include "triangle.hpp"
#include <cmath>

inline double segment_dist_2d(double x, double y, double x0, double y0, double
    x1, double y1)
{
    const double deltaX0 = x - x0;
    const double deltaY0 = y - y0;
    const double deltaX01 = x1 - x0;
    const double deltaY01 = y1 - y0;
    const double lenSq = deltaX01*deltaX01 + deltaY01*deltaY01;
    const double t = std::max(0., std::min(1., (deltaX0*deltaX01 + deltaY0*
        deltaY01)/lenSq));
    const double nearestX = x0 + t*deltaX01;
    const double nearestY = y0 + t*deltaY01;
    return std::sqrt((x - nearestX)*(x - nearestX) + (y - nearestY)*(y -
        nearestY));
}

inline std::array<double, 3> cross(double x0, double y0, double z0, double x1,
    double y1, double z1)
{
    return {y0*z1 - z0*y1, z0*x1 - x0*z1, x0*y1 - y0*x1};
}

inline std::array<double, 3> normalize(const std::array<double, 3>& v)
{
    const double norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return {v[0]/norm, v[1]/norm, v[2]/norm};
}

paz::Triangle::Triangle(double x0, double y0, double z0, double x1, double y1,
    double z1, double x2, double y2, double z2) : x0(x0), y0(y0), z0(z0)
{
    x1 -= x0;
    y1 -= y0;
    z1 -= z0;

    x2 -= x0;
    y2 -= y0;
    z2 -= z0;

    // Guarantees `x1t = y1t = x2t = 0`, `z1t > 0`, `y2t > 0`.
    basisZ = normalize({x1, y1, z1});
    basisX = normalize(cross(x2, y2, z2, basisZ[0], basisZ[1], basisZ[2]));
    basisY = normalize(cross(basisZ[0], basisZ[1], basisZ[2], basisX[0], basisX[1], basisX[2]));

    z1t = basisZ[0]*x1 + basisZ[1]*y1 + basisZ[2]*z1;

    y2t = basisY[0]*x2 + basisY[1]*y2 + basisY[2]*z2;
    z2t = basisZ[0]*x2 + basisZ[1]*y2 + basisZ[2]*z2;
}

// Returns distance from the triangle.
double paz::Triangle::dist(double x, double y, double z) const
{
// SHOULD CHECK FOR COINCIDENT VERTS 1ST AND RETURN APPROPRIATE DIST
    x -= x0;
    y -= y0;
    z -= z0;
    const double xt = basisX[0]*x + basisX[1]*y + basisX[2]*z;
    const double yt = basisY[0]*x + basisY[1]*y + basisY[2]*z;
    const double zt = basisZ[0]*x + basisZ[1]*y + basisZ[2]*z;
    double distInPlane = 0.;
    if(yt <= 0.)
    {
        if(zt <= 0.)
        {
            distInPlane = std::sqrt(yt*yt + zt*zt); // d0
        }
        else if(zt >= z1t)
        {
            const double delta = zt - z1t;
            distInPlane = std::sqrt(yt*yt + delta*delta); // d1
        }
        else
        {
            distInPlane = segment_dist_2d(yt, zt, 0., 0., 0., z1t); // d01
        }
    }
    else
    {
// this is wrong, see collide() below
        const double yRel = yt - y2t;
        const double zRel = zt - z2t;
        if(yRel*-z2t - zRel*-y2t > 0. || yRel*(z1t - z2t) - zRel*-y2t > 0.)
        {
            const double d12 = segment_dist_2d(yt, zt, 0., z1t, y2t, z2t);
            const double d20 = segment_dist_2d(yt, zt, y2t, z2t, 0., 0.);
            distInPlane = std::min(d12, d20);
        }
    }
    return std::sqrt(distInPlane*distInPlane + xt*xt);
}

// Only handles collisions with interior of the triangle (for now).
void paz::Triangle::collide(double x, double y, double z, double& nx, double&
    ny, double& nz, double& d) const
{
    d = 0.;
// SHOULD CHECK FOR COINCIDENT VERTS 1ST AND RETURN EARLY
    x -= x0;
    y -= y0;
    z -= z0;
    const double xt = basisX[0]*x + basisX[1]*y + basisX[2]*z;
    const double yt = basisY[0]*x + basisY[1]*y + basisY[2]*z;
    const double zt = basisZ[0]*x + basisZ[1]*y + basisZ[2]*z;
    if(yt < 0.)
    {
        return;
    }
    if(yt*z2t - zt*y2t > 0)
    {
        return;
    }
    if((zt - z1t)*y2t - yt*(z2t - z1t) > 0) // WRONG
    {
        return;
    }
    /*if(std::abs(xt) >= radius)
    {
        return;
    }*/
    nx = basisX[0];
    ny = basisX[1];
    nz = basisX[2];
    d = xt;
}
