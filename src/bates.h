#include <cmath>
#include <functional>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// From boost
#include <boost/math/quadrature/gauss.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// From alglib
#include "stdafx.h"
#include "interpolation.h"

// Local
#include "fourier.h"

std::complex<long double> CharacteristicFunctionBates(
	std::complex<long double> u, 
	long double tau, long double vol, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda);
	
template <const int precision>
long double EuropeanOptionBates(
	long double tau, long double spot, long double vol, long double strike, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda,
    long double r, long double q);
	
boost::numeric::ublas::vector<double> CalibrateBates(void);
