#ifndef PAZ_ENGINE_UI_HPP
#define PAZ_ENGINE_UI_HPP

#include <PAZ_Engine>
#include <functional>

namespace paz
{
    class Font
    {
        Texture _tex;
        int _scale;
        int _charWidth;

    public:
        Font(Texture tex, int scale, int charWidth);
        Texture tex() const;
        int curScale() const;
        int charWidth() const;
    };

    class Menu;
    class Button
    {
        std::function<std::string(void)> _label;
        std::function<void(Menu&)> _action;
        std::function<bool(void)> _enabled;

    public:
        Button(const std::string& label, const std::function<void(Menu&)>&
            action);
        Button(const std::function<std::string(void)>& label, const std::
            function<void(Menu&)>& action);
        Button(const std::function<std::string(void)>& label, const std::
            function<void(Menu&)>& action, const std::function<bool(void)>&
            enabled);
        std::string label() const;
        void operator()(Menu& m);
        bool enabled() const;
    };

    class Menu
    {
        InstanceBuffer _chars;
        int _curPage;
        int _curButton;

        Font _font;
        std::string _title;
        std::vector<std::vector<Button>> _buttons;

    public:
        Menu(const Font& font, const std::string& title, const std::vector<std::
            vector<Button>>& buttons);
        void update();
        void setState(int page, int button);
        const Font& font() const;
        InstanceBuffer chars() const;
        int curPage() const;
        int curButton() const;
        bool curButtonEnabled() const;
    };
}

#endif
