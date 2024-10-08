#include "PAZ_Engine"
#include "triangle.hpp"
#include "shared.hpp"
#include <limits>
#include <cmath>
#include <numeric>

paz::CollisionMesh::CollisionMesh(const std::string& path, int idx, double
    zOffset, double scale)
{
    std::vector<std::vector<float>> positions;
    std::vector<std::vector<float>> uvs;
    std::vector<std::vector<float>> normals;
    std::vector<std::vector<unsigned int>> indices;
    if(path.size() > 4 && path.substr(path.size() - 4) == ".obj")
    {
        std::vector<std::string> names;
        std::vector<std::vector<unsigned int>> materials;
        std::vector<std::string> materialNames;
        std::vector<std::string> materialLibs;
        parse_obj(get_asset(path), names, positions, uvs, normals, materials,
            materialNames, materialLibs, indices);
    }
    else
    {
        positions.emplace_back();
        uvs.emplace_back();
        normals.emplace_back();
        std::vector<unsigned int> materials;
        std::vector<std::string> materialNames;
        std::vector<std::string> materialLibs;
        indices.emplace_back();
        parse_model(get_asset(path), positions.back(), uvs.back(), normals.
            back(), materials, materialNames, materialLibs, indices.back());
    }
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
}

paz::CollisionMesh::CollisionMesh(const std::vector<float>& positions)
{
    _t = std::make_shared<std::vector<Triangle>>();
    double radiusSq = 0.;
    if(!positions.empty())
    {
        const std::size_t numVertices = positions.size()/4;
        const std::size_t numFaces = numVertices/3;
        _t->reserve(numFaces);
        std::vector<float> normals(4*numVertices, 0.f);
        for(std::size_t i = 0; i < numVertices; i += 3)
        {
            const std::size_t i0 = 4*i;
            const std::size_t i1 = 4*(i + 1);
            const std::size_t i2 = 4*(i + 2);
            const double t0x = positions[i0];
            const double t0y = positions[i0 + 1];
            const double t0z = positions[i0 + 2];
            const double t1x = positions[i1];
            const double t1y = positions[i1 + 1];
            const double t1z = positions[i1 + 2];
            const double t2x = positions[i2];
            const double t2y = positions[i2 + 1];
            const double t2z = positions[i2 + 2];
            _t->emplace_back(t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z);
            std::array<double, 3> nor;
            _t->back().getNormal(nor[0], nor[1], nor[2]);
            for(int j = 0; j < 3; ++j)
            {
                std::copy(nor.begin(), nor.end(), normals.begin() + 4*(i + j));
            }
            radiusSq = std::max(radiusSq, t0x*t0x + t0y*t0y + t0z*t0z);
            radiusSq = std::max(radiusSq, t1x*t1x + t1y*t1y + t1z*t1z);
            radiusSq = std::max(radiusSq, t2x*t2x + t2y*t2y + t2z*t2z);
        }
    }
    _radius = std::sqrt(radiusSq);
}

double paz::CollisionMesh::collide(double x, double y, double z, double radius, double&
    xNew, double& yNew, double& zNew, double& xNor, double& yNor, double& zNor,
    const std::vector<std::size_t>& tris) const
{
    double minDist = inf();
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
    for(auto n : tris)
    {
        double xNorTemp, yNorTemp, zNorTemp, d;
        (*_t)[n].collide(x, y, z, radius, xNorTemp, yNorTemp, zNorTemp, d);
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

void paz::CollisionMesh::castRay(double x, double y, double z, double xDir, double yDir,
    double zDir, double& xNor, double& yNor, double& zNor, double& dist) const
{
    dist = inf();
    xNor = 0.;
    yNor = 0.;
    zNor = 1.;
    // The bounding sphere does not touch the ray.
    const double t = std::max(0., -xDir*x - yDir*y - zDir*z);
    const double deltaX = x + t*xDir;
    const double deltaY = y + t*yDir;
    const double deltaZ = z + t*zDir;
    if(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ > _radius*_radius)
    {
        return;
    }
    // Check each triangle.
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

std::vector<std::size_t> paz::CollisionMesh::sweepVol(const Vec& relPosPrev,
    const Vec& relPos, double radius) const
{
    //TEMP - capped cylinders would give a tighter bound
    const Vec mean = 0.5*(relPosPrev + relPos);
    const double delta = (relPos - mean).norm();
    {
        const double rTotal = delta + _radius + radius;
        const double distSq = relPos.normSq();
        if(distSq > rTotal*rTotal)
        {
            return {};
        }
    }
    std::vector<std::size_t> tris;
    for(std::size_t i = 0; i < _t->size(); ++i)
    {
        Vec mean1(3);
        (*_t)[i].getCentroid(mean1(0), mean1(1), mean1(2));
        const double rTotal = delta + (*_t)[i].radius() + radius;
        const double distSq = (mean1 - mean).normSq();
        if(distSq < rTotal*rTotal)
        {
            tris.push_back(i);
        }
    }
    return tris;
}
