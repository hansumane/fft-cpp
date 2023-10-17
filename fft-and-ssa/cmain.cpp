
#include <iostream>

#include "dft.hpp"
#include "fft.hpp"
#include "ssa.hpp"

using std::cin;
using std::cout;
using std::endl;
using std::string;

int main()
{
    string a, b;

    cout << "a: ";
    cin >> a;
    cout << "b: ";
    cin >> b;

    cout << "a * b = " << ssa::mul(a, b) << endl;

    return 0;
}
