
#include <cmath>

#include "fft.hpp"

using std::abs;
using std::complex;
using std::exp;
using std::round;
using std::size_t;
using std::vector;

vector<complex<double>> get_half(const vector<complex<double>> &arr, bool even)
{
    vector<complex<double>> res;

    for (size_t i = even ? 0 : 1; i < arr.size(); i += 2)
        res.push_back(arr[i]);

    return res;
}

vector<complex<double>> fourier::fft2(const vector<complex<double>> &x)
{
    size_t N = x.size();
    if (N == 1)
        return x;

    vector<complex<double>> xe = fourier::fft2(get_half(x, true));
    vector<complex<double>> xo = fourier::fft2(get_half(x, false));
    vector<complex<double>> factor(N);

    for (size_t i = 0; i < N; i++) {
        complex<double> R = exp(complex<double>(0.0, -2.0) * M_PI * double(i) / double(N));
        factor[i] = xe[i % (N / 2)] + R * xo[i % (N / 2)];
    }

    return factor;
}

vector<complex<double>> _ifft2(const vector<complex<double>> &X)
{
    size_t i, N = X.size();

    if (N == 1)
        return X;

    vector<complex<double>> xe = _ifft2(get_half(X, true));
    vector<complex<double>> xo = _ifft2(get_half(X, false));
    vector<complex<double>> factor(N);

    for (i = 0; i < N; i++) {
        complex<double> R = exp(complex<double>(0.0, 2.0) * M_PI * double(i) / double(N));
        factor[i] = xe[i % (N / 2)] + R * xo[i % (N / 2)];
    }

    return factor;
}

vector<complex<double>> fourier::ifft2(const vector<complex<double>> &X)
{
    vector<complex<double>> x = _ifft2(X);

    for (size_t i = 0; i < x.size(); i++)
        x[i] /= double(x.size());

    return x;
}

vector<int> fourier::ira(const vector<complex<double>> &arr)
{
    vector<int> res(arr.size());

    for (size_t i = 0; i < arr.size(); i++)
        res[i] = static_cast<int>(round(abs(arr[i])));

    return res;
}
