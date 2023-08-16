#include "PAZ_Engine"
#include <cmath>

static constexpr double Speed = 3.;

class Player : public paz::Object
{
public:
    Player() {}
    ~Player() {}
    void update()
    {
        if(paz::Window::KeyDown(paz::Key::W))
        {
            x() -= std::sin(paz::App::CameraYaw())*Speed*paz::Window::FrameTime();
            y() += std::cos(paz::App::CameraYaw())*Speed*paz::Window::FrameTime();
        }
        if(paz::Window::KeyDown(paz::Key::S))
        {
            x() += std::sin(paz::App::CameraYaw())*Speed*paz::Window::FrameTime();
            y() -= std::cos(paz::App::CameraYaw())*Speed*paz::Window::FrameTime();
        }
    }
};

int main(int, char** argv)
{
    const std::string appDir = paz::split_path(argv[0])[0];
    paz::App::Init(appDir + "/scene.frag", {appDir + "/world.obj"});
    Player p;
    paz::App::AttachCamera(p);
    paz::App::Run();
}
