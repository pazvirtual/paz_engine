#include "PAZ_Engine"
#include "triangle.hpp"
#include "shared.hpp"
#include <limits>
#include <cmath>
#include <numeric>

paz::Model::Model(const std::string& path, int idx, double zOffset, double
    scale, const std::string& diffTexPath, const std::array<float, 3>& emiss,
    const std::vector<std::array<double, 9>>& transp) : _emiss(emiss)
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
    _v.addAttribute(4, positions[idx]);
    _v.addAttribute(4, normals[idx]);
    _v.addAttribute(1, std::vector<unsigned int>(positions[idx].size()/4, 1)); //TEMP
    _v.addAttribute(2, uvs[idx]);
    _i = IndexBuffer(indices[idx]);
    if(!diffTexPath.empty())
    {
        _diffTex = Texture(get_asset_image(diffTexPath), MinMagFilter::Linear,
            MinMagFilter::Linear, MipmapFilter::Anisotropic, WrapMode::Repeat,
            WrapMode::Repeat);
    }
    if(!transp.empty())
    {
        std::vector<float> transpPos, transpNor;
        transpPos.reserve(12*transp.size());
        transpNor.reserve(12*transp.size());
        _t->reserve(_t->size() + 2.*transp.size());
        for(const auto& n : transp)
        {
            // Store front and back face because triangles are directional.
            _t->emplace_back(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[
                8]);
            _t->emplace_back(n[6], n[7], n[8], n[3], n[4], n[5], n[0], n[1], n[
                2]);
            radiusSq = std::max(radiusSq, n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            radiusSq = std::max(radiusSq, n[3]*n[3] + n[4]*n[4] + n[5]*n[5]);
            radiusSq = std::max(radiusSq, n[6]*n[6] + n[7]*n[7] + n[8]*n[8]);
            std::array<double, 3> nor;
            _t->back().getNormal(nor[0], nor[1], nor[2]);
            for(int i = 0; i < 3; ++i)
            {
                transpNor.insert(transpNor.end(), nor.begin(), nor.end());
                transpNor.push_back(0.);
                transpPos.insert(transpPos.end(), n.begin() + 3*i, n.begin() +
                    3*i + 3);
                transpPos.push_back(1.);
            }
        }
        _transp.addAttribute(4, transpPos);
        _transp.addAttribute(4, transpNor);
    }
    _radius = std::sqrt(radiusSq);
}

paz::Model::Model(const std::vector<float>& positions, const std::vector<float>&
    uvs, const std::string& diffTexPath, const std::array<float, 3>& emiss,
    const std::vector<std::array<double, 9>>& transp) : _emiss(emiss)
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
        _v.addAttribute(4, positions);
        _v.addAttribute(4, normals);
        _v.addAttribute(1, std::vector<unsigned int>(numVertices, 1));
        _v.addAttribute(2, uvs);
        std::vector<unsigned int> indices(numVertices);
        std::iota(indices.begin(), indices.end(), 0);
        _i = IndexBuffer(indices);
    }
    if(!diffTexPath.empty())
    {
        _diffTex = Texture(get_asset_image(diffTexPath), MinMagFilter::Linear,
            MinMagFilter::Linear, MipmapFilter::Anisotropic, WrapMode::Repeat,
            WrapMode::Repeat);
    }
    if(!transp.empty())
    {
        std::vector<float> transpPos, transpNor;
        transpPos.reserve(12*transp.size());
        transpNor.reserve(12*transp.size());
        _t->reserve(_t->size() + 2.*transp.size());
        for(const auto& n : transp)
        {
            // Store front and back face because triangles are directional.
            _t->emplace_back(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[
                8]);
            _t->emplace_back(n[6], n[7], n[8], n[3], n[4], n[5], n[0], n[1], n[
                2]);
            radiusSq = std::max(radiusSq, n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            radiusSq = std::max(radiusSq, n[3]*n[3] + n[4]*n[4] + n[5]*n[5]);
            radiusSq = std::max(radiusSq, n[6]*n[6] + n[7]*n[7] + n[8]*n[8]);
            std::array<double, 3> nor;
            _t->back().getNormal(nor[0], nor[1], nor[2]);
            for(int i = 0; i < 3; ++i)
            {
                transpNor.insert(transpNor.end(), nor.begin(), nor.end());
                transpNor.push_back(0.);
                transpPos.insert(transpPos.end(), n.begin() + 3*i, n.begin() +
                    3*i + 3);
                transpPos.push_back(1.);
            }
        }
        _transp.addAttribute(4, transpPos);
        _transp.addAttribute(4, transpNor);
    }
    _radius = std::sqrt(radiusSq);
}
