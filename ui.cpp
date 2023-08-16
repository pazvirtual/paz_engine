#include "PAZ_Engine"

void paz::UiDescriptor::setTitle(const std::string& title)
{
    _title = title;
}

void paz::UiDescriptor::setLayout(UiLayout layout)
{
    _layout = layout;
}

void paz::UiDescriptor::alignText(UiAlignment alignment)
{
    _alignment = alignment;
}

void paz::UiDescriptor::addButton(UiAction action, const std::string& str)
{
    _buttons.emplace_back(action, str);
}
