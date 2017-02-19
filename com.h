#include <iostream>
#include <complex>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;
int com()
{
        complex<double> a(2,3), b(4,7);

        Matrix2cd A;

        A<<a,b,
        b,a;

        std::cout<<A<<endl;
}