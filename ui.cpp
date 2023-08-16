#include "ui.hpp"
#include "shared.hpp"

static constexpr double StickSpeed = 10.;

paz::Font::Font(Texture tex, int scale, int charWidth) : _tex(tex), _scale(
    scale), _charWidth(charWidth) {}

paz::Texture paz::Font::tex() const
{
    return _tex;
}

int paz::Font::curScale() const
{
    return std::max(1.f, std::round(_scale*Window::UiScale()));
}

int paz::Font::charWidth() const
{
    return _charWidth;
}

paz::Button::Button(const std::string& label, const std::function<void(Menu&)>&
    action) : _label({[label](){ return label; }}), _action(action), _enabled(
    {[](){ return true; }}) {}

paz::Button::Button(const std::function<std::string(void)>& label, const std::
    function<void(Menu&)>& action) : _label(label), _action(action), _enabled(
    {[](){ return true; }}) {}

paz::Button::Button(const std::function<std::string(void)>& label, const std::
    function<void(Menu&)>& action, const std::function<bool(void)>& enabled) :
    _label(label), _action(action), _enabled(enabled) {}

std::string paz::Button::label() const
{
    return _label();
}

void paz::Button::operator()(Menu& m)
{
    _action(m);
}

bool paz::Button::enabled() const
{
    return _enabled();
}

paz::Menu::Menu(const Font& font, const std::string& title, const std::vector<
    std::vector<Button>>& buttons) : _curPage(0), _curButton(0), _font(font),
    _title(title), _buttons(buttons), _downDist(0.), _upDist(0.) {}

void paz::Menu::update()
{
    if(Window::MouseActive())
    {
        Window::SetCursorMode(CursorMode::Hidden);
    }
    else
    {
        Window::SetCursorMode(CursorMode::Disable);
    }

    if(_curButton >= 0 && (Window::MousePressed(0) || Window::KeyPressed(
        Key::Space) || Window::KeyPressed(Key::Enter) || Window::
        KeyPressed(Key::KeypadEnter) || Window::GamepadPressed(
        GamepadButton::A)))
    {
        _buttons[_curPage][_curButton](*this);
    }
    bool down = Window::KeyPressed(Key::S) || Window::KeyPressed(Key::Down) ||
        Window::GamepadPressed(GamepadButton::Down);
    if(_downDist > 1.)
    {
        down = true;
        _downDist = 0.;
    }
    bool up = Window::KeyPressed(Key::W) || Window::KeyPressed(Key::Up) ||
        Window::GamepadPressed(GamepadButton::Up);
    if(_upDist > 1.)
    {
        up = true;
        _upDist = 0.;
    }
    if(down && !up)
    {
        while(true)
        {
            _curButton = std::min(static_cast<std::size_t>(_curButton) + 1,
                _buttons[_curPage].size() - 1);
            if(static_cast<std::size_t>(_curButton) == _buttons[_curPage].size()
                - 1 || _buttons[_curPage][_curButton].enabled())
            {
                break;
            }
        }
    }
    if(up && !down)
    {
        while(true)
        {
            _curButton = std::max(_curButton - 1, 0);
            if(!_curButton || _buttons[_curPage][_curButton].enabled())
            {
                break;
            }
        }
    }
    const int scale = _font.curScale();
    int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*_font.tex().
        height()));
    const int startRow = (maxVisRows + _buttons[_curPage].size() - 1)/2;

    if(Window::MouseActive())
    {
        _curButton = -1;
        const double x = Window::MousePos().first*Window::DpiScale();
        const double y = Window::MousePos().second*Window::DpiScale();
        for(std::size_t i = 0; i < _buttons[_curPage].size(); ++i)
        {
            const double x0 = 0.;
            const double x1 = x0 + (_buttons[_curPage][i].label().size() + 2)*
                scale*(_font.charWidth() + 1);
            const double y0 = (startRow - i - 1)*scale*_font.tex().height();
            const double y1 = y0 + scale*_font.tex().height();
            if(x >= x0 && x < x1 && y >= y0 && y < y1)
            {
                _curButton = i;
            }
        }
    }

    std::vector<float> highlightAttr;
    std::vector<int> characterAttr;
    std::vector<int> colAttr;
    std::vector<int> rowAttr;
    int row = 0;
    int col = 0;
    for(std::size_t i = 0; i < _buttons[_curPage].size() + 1; ++i)
    {
        float highlight = 0.f;
        if(i && !_buttons[_curPage][i - 1].enabled())
        {
            highlight = -1.f;
        }
        else if(curButtonEnabled() && i == static_cast<std::size_t>(_curButton)
            + 1)
        {
            highlight = 1.f;
        }
        std::string str = (i ? _buttons[_curPage][i - 1].label() : (_curPage ?
            "Options" : _title));
        if(i)
        {
            if(highlight > 0.f)
            {
                str = "->" + str;
            }
            else if(highlight < 0.f)
            {
                str = "  " + str;
            }
            else
            {
                str = "> " + str;
            }
        }
        for(auto n : str)
        {
            if(n == ' ')
            {
                ++col;
            }
            else if(n == '\t')
            {
                col = (col/4 + 1)*4;
            }
            else if(n >= '!' && n <= '~')
            {
                highlightAttr.push_back(highlight);
                characterAttr.push_back(n - '!');
                colAttr.push_back(col);
                rowAttr.push_back(startRow - row);
                ++col;
            }
        }
        ++row;
        col = 0;
    }

    _chars = InstanceBuffer(); //TEMP
    _chars.addAttribute(1, highlightAttr);
    _chars.addAttribute(1, characterAttr);
    _chars.addAttribute(1, colAttr);
    _chars.addAttribute(1, rowAttr);

    const double stickDown = std::max(Window::GamepadLeftStick().second,
        Window::GamepadRightStick().second);
    const double stickUp = std::min(Window::GamepadLeftStick().second, Window::
        GamepadRightStick().second);
    if(stickDown > 0. && stickUp >= 0.)
    {
        _downDist += stickDown*StickSpeed*Window::FrameTime();
    }
    else if(stickUp < 0. && stickDown <= 0.)
    {
        _upDist -= stickUp*StickSpeed*Window::FrameTime();
    }
    else
    {
        _downDist = 0.;
        _upDist = 0.;
    }
}

void paz::Menu::setState(int page, int button)
{
    _curPage = page;
    _curButton = button;
}

const paz::Font& paz::Menu::font() const
{
    return _font;
}

paz::InstanceBuffer paz::Menu::chars() const
{
    return _chars;
}

int paz::Menu::curPage() const
{
    return _curPage;
}

int paz::Menu::curButton() const
{
    return _curButton;
}

bool paz::Menu::curButtonEnabled() const
{
    return _curButton >= 0 && _buttons[_curPage][_curButton].enabled();
}
