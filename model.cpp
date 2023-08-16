#include "PAZ_Engine"
#include "triangle.hpp"
#include <limits>
#include <cmath>

paz::Model::Model(const std::string& path, int idx)
{
    std::vector<std::string> names;
    std::vector<std::vector<float>> positions;
    std::vector<std::vector<float>> uvs;
    std::vector<std::vector<float>> normals;
    std::vector<std::vector<unsigned int>> materials;
    std::vector<std::string> materialNames;
    std::vector<std::string> materialLibs;
    std::vector<std::vector<unsigned int>> indices;
    parse_obj(load_file(path), names, positions, uvs, normals, materials,
        materialNames, materialLibs, indices);
    _t = std::make_shared<std::vector<Triangle>>();
    _t->reserve(indices[idx].size()/3);
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
    }
    _v.attribute(4, positions[idx]);
    _v.attribute(4, normals[idx]);
_v.attribute(1, std::vector<unsigned int>(positions[idx].size()/4, 1));
    _i = IndexBuffer(indices[idx]);
}

double paz::Model::collide(double x, double y, double z, double& gx, double& gy,
    double& gz, double radius, double xPrev, double yPrev, double zPrev, double&
    nx, double& ny, double& nz) const
{
    double minDist = std::numeric_limits<double>::infinity();
    gx = 0.;
    gy = 0.;
    gz = 0.;
    nx = 0.;
    ny = 0.;
    nz = 0.;
const double norm = std::sqrt(x*x + y*y + z*z);
const double rx = x/norm;
const double ry = y/norm;
const double rz = z/norm;
    for(const auto& n : *_t)
    {
        double nxTemp, nyTemp, nzTemp, d;
        n.collide(x, y, z, radius, xPrev, yPrev, zPrev, nxTemp, nyTemp, nzTemp,
            d);
        if(d < radius)
        {
            const double a = radius - d;
            gx += a*nx;
            gy += a*ny;
            gz += a*nz;
            if(d < minDist)
            {
                minDist = d;
                nx = nxTemp;
                ny = nyTemp;
                nz = nzTemp;
            }
        }
    }
const double offset = gx*rx + gy*ry + gz*rz;
std::cout << offset << " " << norm << " " << minDist << std::endl;
if(offset < 0.)
{
throw std::logic_error("NOPE");
}
    return minDist;
}
