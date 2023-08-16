#include "PAZ_Engine"

void paz::InputData::copyEvents(double timestep, double sensitivity)
{
    _timestep = timestep;
    //TEMP - would be nice to copy these structures directly in PAZ_Graphics
    for(int i = 0; i < NumKeys; ++i)
    {
        _keyDown[i] = Window::KeyDown(static_cast<Key>(i));
        _keyPressed[i] = _keyPressed[i] || Window::KeyPressed(static_cast<Key>(
            i));
        _keyReleased[i] = _keyReleased[i] || Window::KeyReleased(static_cast<
            Key>(i));
    }
    for(int i = 0; i < NumMouseButtons; ++i)
    {
        _mouseDown[i] = Window::MouseDown(i);
        _mousePressed[i] = _mousePressed[i] || Window::MousePressed(i);
        _mouseReleased[i] = _mouseReleased[i] || Window::MouseReleased(i);
    }
    _mousePos.first += sensitivity*Window::MousePos().first;
    _mousePos.second += sensitivity*Window::MousePos().second;
    _scrollOffset.first += Window::ScrollOffset().first;
    _scrollOffset.second += Window::ScrollOffset().second;
    for(int i = 0; i < NumGamepadButtons; ++i)
    {
        _gamepadDown[i] = Window::GamepadDown(static_cast<GamepadButton>(i));
        _gamepadPressed[i] = _gamepadPressed[i] || Window::GamepadPressed(
            static_cast<GamepadButton>(i));
        _gamepadReleased[i] = _gamepadReleased[i] || Window::GamepadReleased(
            static_cast<GamepadButton>(i));
    }
    //TEMP - end
    _gamepadLeftStick = Window::GamepadLeftStick();
    _gamepadRightStick = Window::GamepadRightStick();
    _gamepadLeftTrigger = Window::GamepadLeftTrigger();
    _gamepadRightTrigger = Window::GamepadRightTrigger();
    _gamepadActive = Window::GamepadActive();
    _mouseActive = Window::MouseActive();
}

void paz::InputData::resetEvents()
{
    _mousePos = {}; // Assuming cursor is disabled.
    _scrollOffset = {};
    _keyPressed = {};
    _keyReleased = {};
    _mousePressed = {};
    _mouseReleased = {};
    _gamepadPressed = {};
    _gamepadReleased = {};
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
