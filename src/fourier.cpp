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
	long double r, long double q) {
		
	// Create complex numbers and constants.
    std::complex<long double> i (0.0, 1.0);
	long double k = log((spot*exp(-q*tau))/(strike*exp(-r*tau)));
    long double pi = boost::math::constants::pi<long double>();
    
	// Compute integrand.
    auto integrand  = [&] (long double u) {
        return std::real<long double>(exp(i*u*k)*phi(u-(long double)0.5*i))/
			(pow(u, 2.0) + (long double)0.25);
    };
	
	// Initiate integrater and perform integration.
	boost::math::quadrature::gauss<long double, precision> integrator;
    long double I = integrator.integrate(integrand, (long double) 0.0,
		std::numeric_limits<long double>::infinity());
    
    // Compute price.
    long double price = spot*exp(-q*tau) -  (long double) 1.0/pi * 
		sqrt(spot*strike) * exp(-(r+q)*tau / (long double) 2.0) * I;
	
	return price;
}

long double ImpliedVolatilityJackel(long double price, long double tau, 
	long double spot, long double strike, 
	long double r, long double q) {
	
	// Define constants.
	long double deflator = exp(-r*tau);
	long double F = spot*exp((r-q)*tau);
	long double theta = 1;
	
	return implied_volatility_from_a_transformed_rational_guess(price/deflator, F, strike, tau, theta);
}

long double CumulativeFunction(long double x,
	std::function<std::complex<long double> (std::complex<long double>)> phi) {
	// Constants
	std::complex<long double> i (0.0, 1.0);
    long double pi = boost::math::constants::pi<long double>();
	
	// Create integrand.
    auto integrand  = [&] (long double u) {
        return std::real<long double>((exp(-i * u * x) * phi(u)) / (i*u));
    };
	
	// Compute integral.
	boost::math::quadrature::gauss<long double, 1000> integrator;
    long double I = integrator.integrate(integrand, (long double) 0.0,
		std::numeric_limits<long double>::infinity());
		
	// ...
	return 0.5 - 1.0 / pi * I;
}

double BetaConverter(double x) {
    return 2. - 4. * ( sqrt(pow(x, 2.) + x) - x);	
}

std::pair<double, double> LeeParameters(
	std::function<std::complex<long double> (std::complex<long double>)> phi) {
	// Return values.
	double p, q, k, temp;
	double epsilon = 0.0001;
	double delta = 0.0001;
	
	// Create density function.
	auto F  = [&] (double x) {
        return CumulativeFunction(x, phi);
    };
	
		// Compute q.
	q = 0;
	k = 0;
	do {
		temp = q;
		q = -log(F(-k)) / k;
		k = k + delta;
	} while(abs(q - temp) > epsilon);
	
	// Compute p.
	p = 0;
	k = 0;
	do {
		temp = p;
		p = -1. - log(1. - F(k)) / k;
		k = k + delta;
	} while(abs(p - temp) > epsilon);

	return std::pair<double, double> (BetaConverter(p), BetaConverter(q));
 }

// Instantiate template functions.
template long double EuropeanOptionFourier<5000>(std::function<std::complex<long double> (std::complex<long double>)> phi, 
	long double tau, long double spot, long double strike, 
	long double r, long double q);
	
template long double EuropeanOptionFourier<50000>(std::function<std::complex<long double> (std::complex<long double>)> phi, 
	long double tau, long double spot, long double strike, 
	long double r, long double q);
