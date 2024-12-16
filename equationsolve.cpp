#include <iostream>
#include <cmath>
#include <complex>
#include <utility>

std::pair<std::complex<double>, std::complex<double>> solveQuadratic(double a, double b, double c) {
    if (a == 0) {
        throw std::invalid_argument("The coefficient 'a' cannot be zero in a quadratic equation.");
    }

    double D = b * b - 4 * a * c;
    if (D >= 0) {
        std::complex<double> x1 = (-b + std::sqrt(D)) / (2 * a);
        std::complex<double> x2 = (-b - std::sqrt(D)) / (2 * a);
        return std::make_pair(x1, x2);
    } else {
        std::complex<double> x1 = std::complex<double>(-b, std::sqrt(-D)) / (2.0 * a);
        std::complex<double> x2 = std::complex<double>(-b, -std::sqrt(-D)) / (2.0 * a);
        return std::make_pair(x1, x2);
    }
}