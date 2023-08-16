#ifndef PAZ_ENGINE_MATH_HPP
#define PAZ_ENGINE_MATH_HPP

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

inline double fract(const double n)
{
    return n - std::floor(n);
}

inline double normalize_angle(const double n)
{
    return fract(n/(2.*3.14159))*2.*3.14159;
}

#endif
