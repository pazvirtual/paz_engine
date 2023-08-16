#include "PAZ_Engine"

static paz::VertexBuffer SphereV;
static paz::IndexBuffer SphereI;

class Player : public paz::Object {};

class Ball : public paz::Object
{
public:
    Ball()
    {
        v() = SphereV;
        i() = SphereI;
        vis() = true;
        z() = 1.;
        xVel() = 1.;
    }
};

int main(int, char** argv)
{
    const std::string appDir = paz::split_path(argv[0])[0];
    paz::App::Init(appDir + "/scene.frag", {appDir + "/world.obj"});
    Player p;
////
    {
        std::vector<std::string> names;
        std::vector<std::vector<float>> positions;
        std::vector<std::vector<float>> uvs;
        std::vector<std::vector<float>> normals;
        std::vector<std::vector<unsigned int>> materials;
        std::vector<std::string> materialNames;
        std::vector<std::string> materialLibs;
        std::vector<std::vector<unsigned int>> indices;
        paz::parse_obj(paz::load_file(appDir + "/unitsphere.obj"), names,
            positions, uvs, normals, materials, materialNames, materialLibs,
            indices);
        SphereV.attribute(4, positions[0]);
        SphereV.attribute(4, normals[0]);
        // ...
        SphereI = paz::IndexBuffer(indices[0]);
    }
////
    Ball b;
    paz::App::AttachCamera(p);
    paz::App::Run();
}
