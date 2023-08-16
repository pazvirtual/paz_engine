#include "PAZ_Engine"
#include "PAZ_Physics"

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
    _vc = std::make_shared<std::vector<double>>(positions[0].begin(), positions[
        0].end()); // inefficent - 4 vs 3 components
    _ic = std::make_shared<std::vector<std::size_t>>(indices[0].begin(),
        indices[0].end());
    _vg.attribute(4, positions[0]);
    _vg.attribute(4, normals[0]);
    _ig = paz::IndexBuffer(indices[0]);
}

double paz::Model::collide(double x, double y, double z, double& gx, double& gy,
    double& gz, double radius) const
{
    double minDist = std::numeric_limits<double>::infinity();
    gx = 0.;
    gy = 0.;
    gz = 0.;
    for(std::size_t i = 0; i < _ic->size(); i += 3)
    {
        double nx, ny, nz, d;
{
        std::size_t i0 = 4*(*_ic)[i];
        std::size_t i1 = 4*(*_ic)[i + 1];
        std::size_t i2 = 4*(*_ic)[i + 2];
        double t0x = (*_vc)[i0];
        double t0y = (*_vc)[i0 + 1];
        double t0z = (*_vc)[i0 + 2];
        double t1x = (*_vc)[i1];
        double t1y = (*_vc)[i1 + 1];
        double t1z = (*_vc)[i1 + 2];
        double t2x = (*_vc)[i2];
        double t2y = (*_vc)[i2 + 1];
        double t2z = (*_vc)[i2 + 2];
        paz::Triangle t(t0x, t0y, t0z, t1x, t1y, t1z, t2x, t2y, t2z);
        t.collide(x, y, z, nx, ny, nz, d);
}
        const double absD = std::abs(d);
        //const int a = sgn(d);
        const double a = -(radius + d);
        if(d && absD < radius && absD < minDist)
        {
            minDist = absD;
            gx += a*nx;
            gy += a*ny;
            gz += a*nz;
        }
    }
    return minDist;
}
