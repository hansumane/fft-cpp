#ifndef _DFT_HPP_
#define _DFT_HPP_

#include <complex>
#include <vector>

namespace fourier {

std::vector<std::complex<double>> dft(const std::vector<std::complex<double>> &x);
std::vector<std::complex<double>> idft(const std::vector<std::complex<double>> &X);

}

#endif /* _DFT_HPP_ */
