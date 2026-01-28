#include "roots.hpp"
#include <cmath>
#include <functional>

constexpr double TOL = 1e-6;
constexpr int MAX_ITER = 1000000;

/*
 * Bisection method:
 * Repeatedly halves the interval until the root is found
 */
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    double fa = f(a);
    double fb = f(b);

    // Root must be between a and b
    if (fa * fb > 0)
        return false;

    for (int i = 0; i < MAX_ITER; i++)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (std::abs(fc) < TOL)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    *root = 0.5 * (a + b);
    return true;
}

/*
 * Regula falsi (false position):
 * Similar to bisection but uses a straight-line estimate
 */
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
        return false;

    for (int i = 0; i < MAX_ITER; i++)
    {
        double c = b - fb * (b - a) / (fb - fa);
        double fc = f(c);

        if (std::abs(fc) < TOL)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    *root = b - fb * (b - a) / (fb - fa);
    return true;
}

/*
 * Newton–Raphson method:
 * Uses the derivative to improve the guess
 */
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    double x = c;

    for (int i = 0; i < MAX_ITER; i++)
    {
        double fx = f(x);
        double gx = g(x);

        if (std::abs(fx) < TOL)
        {
            *root = x;
            return true;
        }

        // Avoid divide by zero
        if (gx == 0.0)
            return false;

        x = x - fx / gx;
    }

    *root = x;
    return true;
}

/*
 * Secant method:
 * Similar to Newton but does not use the derivative
 */
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    double x0 = a;
    double x1 = b;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < MAX_ITER; i++)
    {
        if (std::abs(f1 - f0) < 1e-12)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (std::abs(f(x2)) < TOL)
        {
            *root = x2;
            return true; 
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    *root = x1;
    return true;
}
