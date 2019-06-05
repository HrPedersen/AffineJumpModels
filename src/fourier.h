#include <functional>
#include <cmath>
#include <complex>

// From Lets be rational article.
#include "lets_be_rational.h"

// From Boost.
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/normal.hpp>

template <const int precision>
long double EuropeanOptionFourier(std::function<std::complex<long double> (std::complex<long double>)> phi, 
	long double tau, long double spot, long double strike, 
	long double r, long double q);

template <const int precision>
long double ImpliedVolatilityFourier(
	std::function<std::complex<long double> (std::complex<long double>)> phi, 
	long double tau, long double spot, long double strike, 
	long double r, long double q);

long double ImpliedVolatilityJackel(long double price, long double tau, 
	long double spot, long double strike, 
	long double r, long double q); 

std::pair<double, double> LeeParameters(
	std::function<std::complex<long double> (std::complex<long double>)> phi);
