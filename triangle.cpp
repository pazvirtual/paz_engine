#include "triangle.hpp"
#include <cmath>
#include <limits>

static constexpr std::size_t NumSteps = 100'000;

inline double square(double x)
{
    return x*x;
}

inline double segment_dist_sq(double x, double y, double x0, double y0, double
    x1, double y1, double& nearestDeltaX, double& nearestDeltaY)
{
    const double deltaX0 = x - x0;
    const double deltaY0 = y - y0;
    const double deltaX01 = x1 - x0;
    const double deltaY01 = y1 - y0;
    const double lenSq = deltaX01*deltaX01 + deltaY01*deltaY01;
    const double t = std::max(0., std::min(1., (deltaX0*deltaX01 + deltaY0*
        deltaY01)/lenSq));
    nearestDeltaX = x - x0 - t*deltaX01;
    nearestDeltaY = y - y0 - t*deltaY01;
    return square(nearestDeltaX) + square(nearestDeltaY);
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
    double z1, double x2, double y2, double z2) : x0(x0), y0(y0), z0(z0),
    _degenerate((x0 == x1 && y0 == y1 && z0 == z1) || (x0 == x2 && y0 == y2 &&
    z0 == z2) || (x1 == x2 && y1 == y2 && z1 == z2)), _centroid({(x0 + x1 + x2)/
    3., (y0 + y1 + y2)/3., (z0 + z1 + z2)/3.})
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
    basisY = normalize(cross(basisZ[0], basisZ[1], basisZ[2], basisX[0], basisX[
        1], basisX[2]));

    z1t = basisZ[0]*x1 + basisZ[1]*y1 + basisZ[2]*z1;

    y2t = basisY[0]*x2 + basisY[1]*y2 + basisY[2]*z2;
    z2t = basisZ[0]*x2 + basisZ[1]*y2 + basisZ[2]*z2;

    _radius = square(x0 - _centroid[0]) + square(y0 - _centroid[1]) + square(z0
        - _centroid[2]);
    _radius = std::min(_radius, square(x1 - _centroid[0]) + square(y1 -
        _centroid[1]) + square(z1 - _centroid[2]));
    _radius = std::min(_radius, square(x2 - _centroid[0]) + square(y2 -
        _centroid[1]) + square(z2 - _centroid[2]));
    _radius = std::sqrt(_radius);
}

double paz::Triangle::dist_transformed(double xt, double yt, double zt, double&
    deltaYt, double& deltaZt) const
{
    // The triangle is degenerate.
    if(_degenerate)
    {
        return std::numeric_limits<double>::infinity(); //TEMP
    }

    const double absXt = std::abs(xt);

    // Check bottom edge.
    if(yt < 0.)
    {
        return std::sqrt(segment_dist_sq(yt, zt, 0., 0., 0., z1t, deltaYt,
            deltaZt) + xt*xt);
    }

    // Check interior.
    const double slopeLeft = y2t/z2t;
    const double yLeft = slopeLeft*zt;
    if(z2t == z1t)
    {
        if(zt < z1t && yt < yLeft)
        {
            return absXt;
        }
    }
    const double slopeRight = y2t/(z2t - z1t);
    const double yRight = slopeRight*(zt - z1t);
    if(!z2t)
    {
        if(zt > 0. && yt < yRight)
        {
            return absXt;
        }
    }
    if(slopeLeft < 0.)
    {
        if(yt > yLeft && yt < yRight)
        {
            return absXt;
        }
    }
    if(slopeRight > 0.)
    {
        if(yt < yLeft && yt > yRight)
        {
            return absXt;
        }
    }
    if(yt < yLeft && yt < yRight)
    {
        return absXt;
    }

    // Check remaining edges.
    double dyzSq = segment_dist_sq(yt, zt, 0., 0., y2t, z2t, deltaYt, deltaZt);
    double deltaYtNew;
    double deltaZtNew;
    const double dyzSqNew = segment_dist_sq(yt, zt, 0., z1t, y2t, z2t,
        deltaYtNew, deltaZtNew);
    if(dyzSqNew < dyzSq)
    {
        dyzSq = dyzSqNew;
        deltaYt = deltaYtNew;
        deltaZt = deltaZtNew;
    }
    return std::sqrt(dyzSq + xt*xt);
}

void paz::Triangle::collide(double x, double y, double z, double radius, double
    xPrev, double yPrev, double zPrev, double& nx, double& ny, double& nz,
    double& d) const
{
    d = std::numeric_limits<double>::infinity();
    nx = 0.;
    ny = 0.;
    nz = 0.;

    // The triangle is degenerate.
    if(_degenerate)
    {
        return; //TEMP
    }

    // The sphere is too far from the triangle
    if(square(x - _centroid[0]) + square(y - _centroid[1]) + square(z -
        _centroid[2]) > square(_radius + radius))
    {
        return;
    }

    x -= x0;
    xPrev -= x0;
    y -= y0;
    z -= z0;
    yPrev -= y0;
    zPrev -= z0;
    const double xt = basisX[0]*x + basisX[1]*y + basisX[2]*z;
    const double xPrevt = basisX[0]*xPrev + basisX[1]*yPrev + basisX[2]*zPrev;

    // The sphere has always been behind the triangle.
    if(xPrevt > 0. && xt > 0.)
    {
        return;
    }

    const double yt = basisY[0]*x + basisY[1]*y + basisY[2]*z;
    const double zt = basisZ[0]*x + basisZ[1]*y + basisZ[2]*z;
    const double yPrevt = basisY[0]*xPrev + basisY[1]*yPrev + basisY[2]*zPrev;
    const double zPrevt = basisZ[0]*xPrev + basisZ[1]*yPrev + basisZ[2]*zPrev;

    const double deltaXt = xt - xPrevt;

    // Find closest approach to the triangle.
    double nearestXt = xt;
    double nearestYt = 0.;
    double nearestZt = 0.;
    if(deltaXt)
    {
        const double ySlope = (yt - yPrevt)/deltaXt;
        const double zSlope = (zt - zPrevt)/deltaXt;
        for(std::size_t i = 0; i < NumSteps; ++i) //TEMP - analytical solution?
        {
            const double deltaXTempt = i*deltaXt/(NumSteps - 1);
            const double xTempt = deltaXTempt + xPrevt;
            const double yTempt = ySlope*deltaXTempt + yPrevt;
            const double zTempt = zSlope*deltaXTempt + zPrevt;
            double nearestYTempt = 0.;
            double nearestZTempt = 0.;
            const double dNew = dist_transformed(xTempt, yTempt, zTempt,
                nearestYTempt, nearestZTempt);
            if(dNew < d)
            {
                d = dNew;
                nearestXt = xTempt;
                nearestYt = nearestYTempt;
                nearestZt = nearestZTempt;
            }
        }
    }
    else
    {
        d = dist_transformed(xt, yt, zt, nearestYt, nearestZt);
    }

    // The sphere has never touched the triangle.
    if(d > radius)
    {
        return;
    }

    double norm = std::sqrt(square(xt) + square(nearestYt) + square(nearestZt));
    if(nearestXt > 0.)
    {
        norm = -norm;
    }
    const double dirX = nearestXt/norm;
    const double dirY = nearestYt/norm;
    const double dirZ = nearestZt/norm;
    nx = basisX[0]*dirX + basisY[0]*dirY + basisZ[0]*dirZ;
    ny = basisX[1]*dirX + basisY[1]*dirY + basisZ[1]*dirZ;
    nz = basisX[2]*dirX + basisY[2]*dirY + basisZ[2]*dirZ;
}
