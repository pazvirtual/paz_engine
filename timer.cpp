#include "timer.hpp"
#include <cmath>

paz::Timer::Timer() : _stopped(false), _initial(std::chrono::steady_clock::
    now()) {}

void paz::Timer::start()
{
    _stopped = false;
    _initial = std::chrono::steady_clock::now();
}

void paz::Timer::stop()
{
    _stopped = true;
    _final = std::chrono::steady_clock::now();
}

double paz::Timer::get() const
{
    if(_stopped)
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(_final -
            _initial).count()*1e-6;
    }
    else
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(std::
            chrono::steady_clock::now() - _initial).count()*1e-6;
    }
}

double paz::Timer::getAndRestart()
{
    double delta;
    if(_stopped)
    {
        delta = std::chrono::duration_cast<std::chrono::microseconds>(_final -
            _initial).count()*1e-6;
        start();
    }
    else
    {
        const auto now = std::chrono::steady_clock::now();
        delta = std::chrono::duration_cast<std::chrono::microseconds>(now -
            _initial).count()*1e-6;
        _initial = now;
    }
    return delta;
}
