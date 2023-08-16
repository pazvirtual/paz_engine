#include "PAZ_Engine"
#include "triangle.hpp"
#include "shared.hpp"
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
    parse_obj(getAsset(path).str(), names, positions, uvs, normals, materials,
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
//TEMP - overwrite UVs b/c sphere50.obj's aren't right
for(std::size_t i = 0; i < positions[idx].size()/4; ++i)
{
const double x = positions[idx][4*i];
const double y = positions[idx][4*i + 1];
const double z = positions[idx][4*i + 2];
const double rad = std::sqrt(x*x + y*y + z*z);
const double lat = std::asin(z/rad);
const double lon = std::atan2(y, x);
uvs[idx][2*i] = 0.5*lon/M_PI + 0.5;
uvs[idx][2*i + 1] = lat/M_PI + 0.5;
}
//TEMP
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
