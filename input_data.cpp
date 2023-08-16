#include "PAZ_Engine"

paz::InputData::InputData(double timestep) : _timestep(timestep)
{
    //TEMP - would be nice to copy these structures directly in PAZ_Graphics
    for(int i = 0; i < NumKeys; ++i)
    {
        _keyDown[i] = Window::KeyDown(static_cast<paz::Key>(i));
        _keyPressed[i] = Window::KeyPressed(static_cast<paz::Key>(i));
        _keyReleased[i] = Window::KeyReleased(static_cast<paz::Key>(i));
    }
    for(int i = 0; i < NumMouseButtons; ++i)
    {
        _mouseDown[i] = Window::MouseDown(i);
        _mousePressed[i] = Window::MousePressed(i);
        _mouseReleased[i] = Window::MouseReleased(i);
    }
    _mousePos = paz::Window::MousePos();
    _scrollOffset = paz::Window::ScrollOffset();
    for(int i = 0; i < NumGamepadButtons; ++i)
    {
        _gamepadDown[i] = Window::GamepadDown(static_cast<paz::GamepadButton>(
            i));
        _gamepadPressed[i] = Window::GamepadPressed(static_cast<paz::
            GamepadButton>(i));
        _gamepadReleased[i] = Window::GamepadReleased(static_cast<paz::
            GamepadButton>(i));
    }
    //TEMP - end
    _gamepadLeftStick = paz::Window::GamepadLeftStick();
    _gamepadRightStick = paz::Window::GamepadRightStick();
    _gamepadLeftTrigger = paz::Window::GamepadLeftTrigger();
    _gamepadRightTrigger = paz::Window::GamepadRightTrigger();
    _gamepadActive = paz::Window::GamepadActive();
    _mouseActive = paz::Window::MouseActive();
}

double paz::InputData::timestep() const
{
    return _timestep;
}

bool paz::InputData::keyDown(Key key) const
{
    return _keyDown.at(static_cast<int>(key));
}

bool paz::InputData::keyPressed(Key key) const
{
    return _keyPressed.at(static_cast<int>(key));
}

bool paz::InputData::keyReleased(Key key) const
{
    return _keyReleased.at(static_cast<int>(key));
}

bool paz::InputData::mouseDown(int button) const
{
    return _mouseDown.at(button);
}

bool paz::InputData::mousePressed(int button) const
{
    return _mousePressed.at(button);
}

bool paz::InputData::mouseReleased(int button) const
{
    return _mouseReleased.at(button);
}

const std::pair<double, double>& paz::InputData::mousePos() const
{
    return _mousePos;
}

const std::pair<double, double>& paz::InputData::scrollOffset() const
{
    return _scrollOffset;
}

bool paz::InputData::gamepadDown(GamepadButton button) const
{
    return _gamepadDown.at(static_cast<int>(button));
}

bool paz::InputData::gamepadPressed(GamepadButton button) const
{
    return _gamepadPressed.at(static_cast<int>(button));
}

bool paz::InputData::gamepadReleased(GamepadButton button) const
{
    return _gamepadReleased.at(static_cast<int>(button));
}

const std::pair<double, double>& paz::InputData::gamepadLeftStick() const
{
    return _gamepadLeftStick;
}

const std::pair<double, double>& paz::InputData::gamepadRightStick() const
{
    return _gamepadRightStick;
}

double paz::InputData::gamepadLeftTrigger() const
{
    return _gamepadLeftTrigger;
}

double paz::InputData::gamepadRightTrigger() const
{
    return _gamepadRightTrigger;
}

bool paz::InputData::gamepadActive() const
{
    return _gamepadActive;
}

bool paz::InputData::mouseActive() const
{
    return _mouseActive;
}
