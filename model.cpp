#include "PAZ_Engine"
#include "triangle.hpp"
#include <limits>

paz::Model::Model() {}

paz::Model::Model(const std::string& path)
{
    std::vector<std::string> names;
    std::vector<std::vector<float>> positions;
    std::vector<std::vector<float>> uvs;
    std::vector<std::vector<float>> normals;
    std::vector<std::vector<unsigned int>> materials;
    std::vector<std::string> materialNames;
    std::vector<std::string> materialLibs;
    std::vector<std::vector<unsigned int>> indices;
    paz::parse_obj(paz::load_file(path), names, positions, uvs, normals,
        materials, materialNames, materialLibs, indices);
    _t = std::make_shared<std::vector<paz::Triangle>>();
    _t->reserve(indices[0].size()/3);
    for(std::size_t i = 0; i < indices[0].size(); i += 3)
    {
        const std::size_t i0 = 4*indices[0][i];
        const std::size_t i1 = 4*indices[0][i + 1];
        const std::size_t i2 = 4*indices[0][i + 2];
        const double t0x = positions[0][i0];
        const double t0y = positions[0][i0 + 1];
        const double t0z = positions[0][i0 + 2];
        const double t1x = positions[0][i1];
        const double t1y = positions[0][i1 + 1];
        const double t1z = positions[0][i1 + 2];
        const double t2x = positions[0][i2];
        const double t2y = positions[0][i2 + 1];
        const double t2z = positions[0][i2 + 2];
        _t->emplace_back(t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z);
    }
    _v.attribute(4, positions[0]);
    _v.attribute(4, normals[0]);
    _i = paz::IndexBuffer(indices[0]);
}

double paz::Model::collide(double x, double y, double z, double& gx, double& gy,
    double& gz, double radius, double xPrev, double yPrev, double zPrev) const
{
    double minDist = std::numeric_limits<double>::infinity();
    gx = 0.;
    gy = 0.;
    gz = 0.;
    for(const auto& n : *_t)
    {
        double nx, ny, nz, d;
        n.collide(x, y, z, radius, xPrev, yPrev, zPrev, nx, ny, nz, d);
        if(d < radius && d < minDist)
        {
            minDist = d;
            gx += nx;
            gy += ny;
            gz += nz;
        }
    }
    return minDist;
}
