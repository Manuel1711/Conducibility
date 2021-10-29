#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include <iostream>

using namespace std;

namespace bm=boost::multiprecision;
namespace bq=boost::math::quadrature;

using Real=
   bm::number<bm::mpfr_float_backend<128>>;

int main()
{
   const Real inf=
     std::numeric_limits<Real>::infinity();

   const auto f=
     [](const Real& x) -> Real
     {
       return exp(-x*x/2);
     };

   const Real infLimit=0.23532;
   const Real supLimit=inf;

   const Real integr=
     bq::gauss_kronrod<Real,61>::integrate(f,infLimit,supLimit,5,1e-16);

   cout<<integr<<endl;

   return 0;
}
