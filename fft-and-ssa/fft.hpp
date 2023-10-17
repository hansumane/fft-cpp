#ifndef _FFT_HPP_
#define _FFT_HPP_

#include <complex>
#include <vector>

namespace fourier {

std::vector<std::complex<double>> fft2(const std::vector<std::complex<double>> &x);
std::vector<std::complex<double>> ifft2(const std::vector<std::complex<double>> &X);
std::vector<int> ira(const std::vector<std::complex<double>> &arr);

}


#endif /* _FFT_HPP_ */
