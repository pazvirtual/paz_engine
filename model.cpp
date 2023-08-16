#include "PAZ_Engine"
#include "triangle.hpp"
#include "shared.hpp"
#include <limits>
#include <cmath>

paz::Model::Model(const std::string& path, int idx, double zOffset, double
    scale)
{
    std::vector<std::string> names;
    std::vector<std::vector<float>> positions;
    std::vector<std::vector<float>> uvs;
    std::vector<std::vector<float>> normals;
    std::vector<std::vector<unsigned int>> materials;
    std::vector<std::string> materialNames;
    std::vector<std::string> materialLibs;
    std::vector<std::vector<unsigned int>> indices;
    parse_obj(get_asset(path).str(), names, positions, uvs, normals, materials,
        materialNames, materialLibs, indices);
    _t = std::make_shared<std::vector<Triangle>>();
    _t->reserve(indices[idx].size()/3);
    if(std::abs(zOffset) > 1e-6)
    {
        for(std::size_t i = 0; i < positions[idx].size(); i += 4)
        {
            positions[idx][i + 2] += zOffset;
        }
    }
    if(scale > 0. && std::abs(scale - 1.) > 1e-6)
    {
        for(std::size_t i = 0; i < positions[idx].size(); i += 4)
        {
            for(std::size_t j = 0; j < 3; ++j)
            {
                positions[idx][i + j] *= scale;
            }
        }
    }
    double radiusSq = 0.;
    for(std::size_t i = 0; i < indices[idx].size(); i += 3)
    {
        const std::size_t i0 = 4*indices[idx][i];
        const std::size_t i1 = 4*indices[idx][i + 1];
        const std::size_t i2 = 4*indices[idx][i + 2];
        const double t0x = positions[idx][i0];
        const double t0y = positions[idx][i0 + 1];
        const double t0z = positions[idx][i0 + 2];
        const double t1x = positions[idx][i1];
        const double t1y = positions[idx][i1 + 1];
        const double t1z = positions[idx][i1 + 2];
        const double t2x = positions[idx][i2];
        const double t2y = positions[idx][i2 + 1];
        const double t2z = positions[idx][i2 + 2];
        _t->emplace_back(t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z);
        radiusSq = std::max(radiusSq, t0x*t0x + t0y*t0y + t0z*t0z);
        radiusSq = std::max(radiusSq, t1x*t1x + t1y*t1y + t1z*t1z);
        radiusSq = std::max(radiusSq, t2x*t2x + t2y*t2y + t2z*t2z);
    }
    _radius = std::sqrt(radiusSq);
    _v.attribute(4, positions[idx]);
    _v.attribute(4, normals[idx]);
    _v.attribute(1, std::vector<unsigned int>(positions[idx].size()/4, 1)); //TEMP
    _v.attribute(2, uvs[idx]);
    _i = IndexBuffer(indices[idx]);
}

double paz::Model::collide(double x, double y, double z, double radius, double&
    xNew, double& yNew, double& zNew, double& xNor, double& yNor, double& zNor)
    const
{
    double minDist = std::numeric_limits<double>::infinity();
    xNor = 0.;
    yNor = 0.;
    zNor = 1.;
    if(std::sqrt(x*x + y*y + z*z) > _radius + radius)
    {
        return minDist;
    }
    double gx = 0.;
    double gy = 0.;
    double gz = 0.;
    for(const auto& n : *_t)
    {
        double xNorTemp, yNorTemp, zNorTemp, d;
        n.collide(x, y, z, radius, xNorTemp, yNorTemp, zNorTemp, d);
        if(d < radius)
        {
            const double a = radius - d;
            gx += a*xNorTemp;
            gy += a*yNorTemp;
            gz += a*zNorTemp;
            if(d < minDist)
            {
                minDist = d;
                xNor = xNorTemp;
                yNor = yNorTemp;
                zNor = zNorTemp;
            }
        }
    }
    xNew = x + gx;
    yNew = y + gy;
    zNew = z + gz;
    return minDist;
}

void paz::Model::castRay(double x, double y, double z, double xDir, double yDir,
    double zDir, double& xNor, double& yNor, double& zNor, double& dist) const
{
    dist = std::numeric_limits<double>::infinity();
    for(const auto& n : *_t)
    {
        const double d = n.castRay(x, y, z, xDir, yDir, zDir);
        if(d < dist)
        {
            dist = d;
            n.getNormal(xNor, yNor, zNor);
        }
    }
}

bool paz::Model::sweepVol(double xPrev0, double yPrev0, double zPrev0, double
    x0, double y0, double z0, double xPrev1, double yPrev1, double zPrev1,
    double x1, double y1, double z1, double radius) const
{
    //TEMP - capped cylinders would give a tighter bound
    const double xMean0 = 0.5*(xPrev0 + x0);
    const double yMean0 = 0.5*(yPrev0 + y0);
    const double zMean0 = 0.5*(zPrev0 + z0);
    const double delta0 = std::sqrt((x0 - xMean0)*(x0 - xMean0) + (y0 - yMean0)*
        (y0 - yMean0) + (z0 - zMean0)*(z0 - zMean0));
    const double xMean1 = 0.5*(xPrev1 + x1);
    const double yMean1 = 0.5*(yPrev1 + y1);
    const double zMean1 = 0.5*(zPrev1 + z1);
    const double delta1 = std::sqrt((x1 - xMean1)*(x1 - xMean1) + (y1 - yMean1)*
        (y1 - yMean1) + (z1 - zMean1)*(z1 - zMean1));
    const double rTotal = delta0 + delta1 + _radius + radius;
    const double distSq = (xMean1 - xMean0)*(xMean1 - xMean0) + (yMean1 -
        yMean0)*(yMean1 - yMean0) + (zMean1 - zMean0)*(zMean1 - zMean0);
    return distSq < rTotal*rTotal;
}
