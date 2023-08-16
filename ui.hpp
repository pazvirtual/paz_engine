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
        friend class Menu;

        Menu* _parent;
        int _mode;
        std::vector<std::string> _labels;
        std::function<void(Button&)> _action;

    public:
        Button(const std::vector<std::string>& labels, const std::function<void(
            Button&)>& action);
        int mode() const;
        void setMode(int mode);
        const std::string& label() const;
        void operator()();
        Menu& parent() const;
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
    };
}

#endif
