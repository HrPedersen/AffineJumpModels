#include <cmath>
#include <functional>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// From boost
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/quadrature/gauss.hpp>

// From alglib
#include "stdafx.h"
#include "interpolation.h"

// Local
#include "fourier.h"
#include "black.h"

// Functions
std::complex<long double> CharacteristicFunctionMerton(
	std::complex<long double> u, long double tau, 
	long double sigma, long double alpha, long double delta, 
	long double lambda);
	
template <const int precision>
long double EuropeanOptionMerton(long double tau, long double spot, 
	long double strike, long double sigma, long double alpha, long double delta,
	long double lambda, long double r, long double q);
	
boost::numeric::ublas::vector<double> CalibrateMerton(void);
