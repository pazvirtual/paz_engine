#ifndef PAZ_ENGINE_TIMER_HPP
#define PAZ_ENGINE_TIMER_HPP

#include <chrono>
#include <iostream>

namespace paz
{
    class Timer
    {
        bool _stopped;
        std::chrono::time_point<std::chrono::steady_clock> _initial;
        std::chrono::time_point<std::chrono::steady_clock> _final;

    public:
        Timer();
        void start();
        void stop();
        double get() const;
        double getAndRestart();
    };
}

#endif
