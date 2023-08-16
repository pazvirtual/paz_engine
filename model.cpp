#include "PAZ_Engine"
#include "triangle.hpp"
#include "shared.hpp"
#include <limits>
#include <cmath>

static constexpr std::size_t NumSteps = 100;

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
    _v.attribute(4, positions[idx]);
    _v.attribute(4, normals[idx]);
    _v.attribute(1, std::vector<unsigned int>(positions[idx].size()/4, 1)); //TEMP
    _i = IndexBuffer(indices[idx]);
}

double paz::Model::collide(double x, double y, double z, double xPrev, double
    yPrev, double zPrev, double radius, double& xNew, double& yNew, double&
    zNew, double& xNor, double& yNor, double& zNor) const
{
    double minDist = std::numeric_limits<double>::infinity();
    xNew = x;
    yNew = y;
    zNew = z;
    xNor = 0.;
    yNor = 0.;
    zNor = 1.;
bool done = false;
for(std::size_t i = 0; i < NumSteps; ++i)
{
double gx = 0.;
double gz = 0.;
double gy = 0.;
    const double curX = xPrev + (i + 1)*(x - xPrev)/NumSteps;
    const double curY = yPrev + (i + 1)*(y - yPrev)/NumSteps;
    const double curZ = zPrev + (i + 1)*(z - zPrev)/NumSteps;
    for(const auto& n : *_t)
    {
        double xNorTemp, yNorTemp, zNorTemp, d;
        n.collide(curX, curY, curZ, radius, xNorTemp, yNorTemp, zNorTemp, d);
        if(d < radius)
        {
            const double a = radius - d;
            gx += a*xNorTemp;
            gy += a*yNorTemp;
            gz += a*zNorTemp;
            if(d < minDist)
            {
done = true;
                minDist = d;
                xNor = xNorTemp;
                yNor = yNorTemp;
                zNor = zNorTemp;
            }
        }
if(done)
{
xNew = curX + gx;
yNew = curY + gy;
zNew = curZ + gz;
return minDist;
}
    }
}
    return minDist;
}
