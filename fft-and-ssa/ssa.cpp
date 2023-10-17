
#include <cmath>

#include "fft.hpp"
#include "ssa.hpp"

using std::ceil;
using std::complex;
using std::log2;
using std::max;
using std::pow;
using std::string;
using std::vector;

vector<complex<double>> atocv(const string &str)
{
    vector<complex<double>> res;

    /* converting string into complex vector:             *
     * "1234" -> {4.0+0.0i, 3.0+0.0i, 2.0+0.0i, 1.0+0.0i} */
    for (long i = str.size() - 1; i >= 0; i--)
        res.push_back(complex<double>(str[i] - '0', 0.0));

    return res;
}

string cvtoa(const vector<complex<double>> &arr)
{
    long i;
    string res;
    vector<int> temp = fourier::ira(arr);

    /* recombination (carrying): *
     * 0 15 12 0 -> 0 5 3 1      *
     * 0 1> 1>                   */
    for (i = 0; i < long(temp.size() - 1); i++) {
        temp[i + 1] += temp[i] / 10;
        temp[i] %= 10;
    }

    /* skipping zeros:  *
     * 0 5 3 1 0 0 0 0  *
     *               ^i *
     * 0 5 3 1 0 0 0 0  *
     *       ^i         */
    for (i = temp.size() - 1; i >= 0; i--)
        if (temp[i] != 0)
            break;

    while (i >= 0)
        res.push_back(temp[i--] + '0');

    return res;
}

string ssa::mul(const string &a, const string &b)
{
    size_t i;
    vector<complex<double>> vca = atocv(a);
    vector<complex<double>> vcb = atocv(b);

    /* calculating required size (rs):               *
     * 123 (3 ^= 4) -> 0123 (4 * 2 = 8) -> 000001234 */
    size_t ms = max(vca.size(), vcb.size());
    size_t rs = size_t(pow(2.0, ceil(log2(ms)) + 1.0));

    while (vca.size() < rs)
        vca.push_back(0);

    while (vcb.size() < rs)
        vcb.push_back(0);

    vector<complex<double>> fa = fourier::fft2(vca);
    vector<complex<double>> fb = fourier::fft2(vcb);
    vector<complex<double>> fr(rs);

    for (i = 0; i < rs; i++)
        fr[i] = fa[i] * fb[i];

    vector<complex<double>> ifr = fourier::ifft2(fr);

    return cvtoa(ifr);
}
