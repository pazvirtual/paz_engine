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

    class Menu
    {
        InstanceBuffer _chars;
        int _curPage;
        int _curButton;

        Font _font;
        std::string _title;
        std::vector<std::vector<std::pair<std::string, std::function<void(
            Menu&)>>>> _buttons;

    public:
        Menu(const Font& font, const std::string& title, const std::vector<std::
            vector<std::pair<std::string, std::function<void(Menu&)>>>>&
            buttons);
        void update();
        void setState(int page, int button);
        const Font& font() const;
        InstanceBuffer chars() const;
        int curPage() const;
    };
}

#endif
