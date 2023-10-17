
#include <cmath>

#include "dft.hpp"

using std::complex;
using std::exp;
using std::size_t;
using std::vector;

vector<complex<double>> fourier::dft(const vector<complex<double>> &x)
{
    double knN;
    size_t k, n, N = x.size();

    complex<double> res;
    vector<complex<double>> X(N);

    for (k = 0; k < N; k++) {
        res = {0.0, 0.0};
        for (n = 0; n < N; n++) {
            knN = double(k) * double(n) / double(N);
            res += x[n] * exp(complex<double>(0.0, -2.0) * M_PI * knN);
        }
        X[k] = res;
    }

    return X;
}

std::vector<std::complex<double>> fourier::idft(const std::vector<std::complex<double>> &X)
{
    double knN;
    size_t k, n, N = X.size();

    complex<double> res;
    vector<complex<double>> x(N);

    for (n = 0; n < N; n++) {
        res = {0.0, 0.0};
        for (k = 0; k < N; k++) {
            knN = double(k) * double(n) / double(N);
            res += X[k] * exp(complex<double>(0.0, 2.0) * M_PI * knN);
        }
        x[n] = res / double(N);
    }

    return x;
}
