#include "triangle.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

constexpr double square(double x)
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

inline void cross(double x0, double y0, double z0, double x1, double y1, double
    z1, double& xc, double& yc, double& zc)
{
    xc = y0*z1 - z0*y1;
    yc = z0*x1 - x0*z1;
    zc = x0*y1 - y0*x1;
}

inline void normalize(double& x, double& y, double& z)
{
    const double invNorm = 1./std::sqrt(x*x + y*y + z*z);
    x *= invNorm;
    y *= invNorm;
    z *= invNorm;
}

inline bool approx(double a, double b)
{
    return std::abs(a - b) < 1e-6;
}

paz::Triangle::Triangle(double x0, double y0, double z0, double x1, double y1,
    double z1, double x2, double y2, double z2) : _x0(x0), _y0(y0), _z0(z0),
    _degenerate((approx(_x0, x1) && approx(_y0, y1) && approx(_z0, z1)) ||
    (approx(_x0, x2) && approx(_y0, y2) && approx(_z0, z2)) || (approx(x1, x2)
    && approx(y1, y2) && approx(z1, z2))), _centroidX((_x0 + x1 + x2)/3.),
    _centroidY((_y0 + y1 + y2)/3.), _centroidZ((_z0 + z1 + z2)/3.)
{
    x1 -= _x0;
    y1 -= _y0;
    z1 -= _z0;

    x2 -= _x0;
    y2 -= _y0;
    z2 -= _z0;

    // Guarantees `x1t = y1t = x2t = 0`, `z1t > 0`, `y2t > 0`.
    _basisZX = x1;
    _basisZY = y1;
    _basisZZ = z1;
    normalize(_basisZX, _basisZY, _basisZZ);
    cross(x2, y2, z2, _basisZX, _basisZY, _basisZZ, _basisXX, _basisXY,
        _basisXZ);
    normalize(_basisXX, _basisXY, _basisXZ);
    cross(_basisZX, _basisZY, _basisZZ, _basisXX, _basisXY, _basisXZ, _basisYX,
        _basisYY, _basisYZ);
    normalize(_basisYX, _basisYY, _basisYZ);

    _z1t = _basisZX*x1 + _basisZY*y1 + _basisZZ*z1;

    _y2t = _basisYX*x2 + _basisYY*y2 + _basisYZ*z2;
    _z2t = _basisZX*x2 + _basisZY*y2 + _basisZZ*z2;

    _radius = square(_x0 - _centroidX) + square(_y0 - _centroidY) + square(_z0 -
        _centroidZ);
    _radius = std::min(_radius, square(x1 - _centroidX) + square(y1 -
        _centroidY) + square(z1 - _centroidZ));
    _radius = std::min(_radius, square(x2 - _centroidX) + square(y2 -
        _centroidY) + square(z2 - _centroidZ));
    _radius = std::sqrt(_radius);
}

double paz::Triangle::distTransformed(double xt, double yt, double zt, double&
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
        return std::sqrt(segment_dist_sq(yt, zt, 0., 0., 0., _z1t, deltaYt,
            deltaZt) + xt*xt);
    }

    // Check interior.
    const double slopeLeft = _y2t/_z2t;
    const double yLeft = slopeLeft*zt;
    if(approx(_z2t, _z1t))
    {
        if(zt < _z1t && yt < yLeft)
        {
            return absXt;
        }
    }
    else
    {
        const double slopeRight = _y2t/(_z2t - _z1t);
        const double yRight = slopeRight*(zt - _z1t);
        if(!_z2t)
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
        if(slopeLeft > 0. && slopeRight < 0.)
        {
            if(yt < yLeft && yt < yRight)
            {
                return absXt;
            }
        }
    }

    // Check remaining edges.
    double dyzSq = segment_dist_sq(yt, zt, 0., 0., _y2t, _z2t, deltaYt,
        deltaZt);
    double deltaYtNew;
    double deltaZtNew;
    const double dyzSqNew = segment_dist_sq(yt, zt, 0., _z1t, _y2t, _z2t,
        deltaYtNew, deltaZtNew);
    if(dyzSqNew < dyzSq)
    {
        dyzSq = dyzSqNew;
        deltaYt = deltaYtNew;
        deltaZt = deltaZtNew;
    }
    return std::sqrt(dyzSq + xt*xt);
}

void paz::Triangle::collide(double x, double y, double z, double radius, double&
    nx, double& ny, double& nz, double& d) const
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
    if(square(x - _centroidX) + square(y - _centroidY) + square(z - _centroidZ)
        > square(_radius + radius))
    {
        return;
    }

    x -= _x0;
    y -= _y0;
    z -= _z0;
    const double xt = _basisXX*x + _basisXY*y + _basisXZ*z;

    // The sphere is behind the triangle.
    if(xt > 0.)
    {
        return;
    }

    const double yt = _basisYX*x + _basisYY*y + _basisYZ*z;
    const double zt = _basisZX*x + _basisZY*y + _basisZZ*z;

    // Find distance from the triangle.
    double nearestXt = xt;
    double nearestYt = 0.;
    double nearestZt = 0.;
    d = distTransformed(xt, yt, zt, nearestYt, nearestZt);

    // The sphere is not touching the triangle.
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
    nx = _basisXX*dirX + _basisYX*dirY + _basisZX*dirZ;
    ny = _basisXY*dirX + _basisYY*dirY + _basisZY*dirZ;
    nz = _basisXZ*dirX + _basisYZ*dirY + _basisZZ*dirZ;
}

double paz::Triangle::castRay(double x, double y, double z, double xDir, double
    yDir, double zDir) const
{
    // The triangle is degenerate.
    if(_degenerate)
    {
        return std::numeric_limits<double>::infinity();
    }

    // Ray is parallel to or points away from triangle.
    const double rayDotNor = -xDir*_basisXX - yDir*_basisXY - zDir*_basisXZ;
    if(rayDotNor >= 0.)
    {
        return std::numeric_limits<double>::infinity();
    }

    // Point is behind triangle.
    x -= _x0;
    y -= _y0;
    z -= _z0;
    const double posDotNor = -x*_basisXX - y*_basisXY - z*_basisXZ;
    if(posDotNor < 0.)
    {
        return std::numeric_limits<double>::infinity();
    }

    // Compute ray-plane intersection.
    const double dist = -posDotNor/rayDotNor;
    const double xp = x + dist*xDir;
    const double yp = y + dist*yDir;
    const double zp = z + dist*zDir;

    // Finally, check if the intersection is inside of the triangle.
    const double ypt = _basisYX*xp + _basisYY*yp + _basisYZ*zp;
    const double zpt = _basisZX*xp + _basisZY*yp + _basisZZ*zp;

    // Check y limits.
    if(ypt < 0. || ypt > _y2t)
    {
        return std::numeric_limits<double>::infinity();
    }

    // Check interior.
    const double slopeLeft = _z2t/_y2t;
    const double zLeft = slopeLeft*ypt;
    const double slopeRight = (_z2t - _z1t)/_y2t;
    const double zRight = _z1t + slopeRight*ypt;
    if(zpt < zLeft || zpt > zRight)
    {
        return std::numeric_limits<double>::infinity();
    }

    return dist;
}

void paz::Triangle::getNormal(double& xNor, double& yNor, double& zNor) const
{
    xNor = -_basisXX;
    yNor = -_basisXY;
    zNor = -_basisXZ;
}

void paz::Triangle::getCentroid(double& x, double& y, double& z) const
{
    x = _centroidX;
    y = _centroidY;
    z = _centroidZ;
}

double paz::Triangle::radius() const
{
    return _radius;
}
