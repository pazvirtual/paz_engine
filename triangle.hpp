#ifndef PAZ_ENGINE_TRIANGLE_HPP
#define PAZ_ENGINE_TRIANGLE_HPP

namespace paz
{
    class Triangle
    {
        const double _x0, _y0, _z0;
        const bool _degenerate;
        double _basisXX, _basisXY, _basisXZ;
        double _basisYX, _basisYY, _basisYZ;
        double _basisZX, _basisZY, _basisZZ;
        double _z1t, _y2t, _z2t;
        double _centroidX, _centroidY, _centroidZ;
        double _radius;

        double distTransformed(double xt, double yt, double zt, double& deltaYt,
            double& deltaZt) const;

    public:
        Triangle(double x0, double y0, double z0, double x1, double y1, double
            z1, double x2, double y2, double z2);
        void collide(double x, double y, double z, double radius, double& nx,
            double& ny, double& nz, double& d) const;
        double castRay(double x, double y, double z, double xDir, double yDir,
            double zDir) const;
        void getNormal(double& x, double& y, double& z) const;
        void getCentroid(double& x, double& y, double& z) const;
        double radius() const;
    };
}

#endif
