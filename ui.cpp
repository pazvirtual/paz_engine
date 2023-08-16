#include "ui.hpp"
#include "shared.hpp"

paz::Font::Font(Texture tex, int scale, int charWidth) : _tex(tex), _scale(
    scale), _charWidth(charWidth) {}

paz::Texture paz::Font::tex() const
{
    return _tex;
}

int paz::Font::curScale() const
{
    return std::round(_scale*Window::UiScale());
}

int paz::Font::charWidth() const
{
    return _charWidth;
}

paz::Button::Button(const std::string& label, const std::function<void(Menu&)>&
    action) : _label({[label](){ return label; }}), _action(action) {}

paz::Button::Button(const std::function<std::string(void)>& label, const std::
    function<void(Menu&)>& action) : _label(label), _action(action) {}

std::string paz::Button::label() const
{
    return _label();
}

void paz::Button::operator()(Menu& m)
{
    _action(m);
}

paz::Menu::Menu(const Font& font, const std::string& title, const std::vector<
    std::vector<Button>>& buttons) : _curPage(0), _curButton(0), _font(font),
    _title(title), _buttons(buttons) {}

void paz::Menu::update()
{
    if(Window::MouseActive())
    {
        Window::SetCursorMode(CursorMode::Normal);
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
    if(Window::KeyPressed(Key::S) || Window::KeyPressed(Key::Down) ||
        Window::GamepadPressed(GamepadButton::Down))
    {
        _curButton = std::min(static_cast<std::size_t>(_curButton) + 1,
            _buttons[_curPage].size() - 1);
    }
    if(Window::KeyPressed(Key::W) || Window::KeyPressed(Key::Up) ||
         Window::GamepadPressed(GamepadButton::Up))
    {
        _curButton = std::max(_curButton - 1, 0);
    }
    const int scale = _font.curScale();
    int maxVisRows = std::ceil(Window::ViewportHeight()/(scale*_font.tex().
        height()));
    const int startRow = (maxVisRows + _buttons[_curPage].size() - 1)/2;

    if(Window::MouseActive())
    {
        _curButton = -1;
        for(std::size_t i = 0; i < _buttons[_curPage].size(); ++i)
        {
            const double x0 = 0;
            const double x1 = x0 + (_buttons[_curPage][i].label().size() + 2)*
                scale*(_font.charWidth() + 1);
            const double y0 = (startRow - i - 1)*scale*_font.tex().height();
            const double y1 = y0 + scale*_font.tex().height();
            if(Window::MousePos().first >= x0 && Window::MousePos().
                first < x1 && Window::MousePos().second >= y0 &&
                Window::MousePos().second < y1)
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
        const bool highlight = _curButton >= 0 && i == static_cast<std::
            size_t>(_curButton) + 1;
        std::string str = (i ? _buttons[_curPage][i - 1].label() : (_curPage ?
            "Options" : _title));
        if(i)
        {
            if(highlight)
            {
                str = "->" + str;
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
