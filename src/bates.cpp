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

std::complex<long double> CharacteristicFunctionBates(
	std::complex<long double> u, 
	long double tau, long double vol, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda) {
	
    // Returns the characteristic function of X_T conditional on x_t under the 
    // Q-measure. 
	
	// Create complex numbers and constants.
    std::complex<long double> i (0.0, 1.0);
	
	// First compute the Heston coefficients.
	// Declare d.
    std::complex<long double> d = sqrt(pow(rho*sigma*u*i - kappa, 2.0) - 
        pow(sigma, 2.0)*(-u*i - pow(u, 2.0)));
        
    // Declare c.
    std::complex<long double> c = (kappa - rho*sigma*i*u - d)/
        (kappa - rho*sigma*i*u + d);
		
	// Declare A.
	std::complex<long double> A = theta*kappa/pow(sigma, 2.0)*
		((kappa - rho*sigma*u*i - d)*tau -
		(long double)2.0*log(((long double)1.0-c*exp(-d*tau))/
		((long double)1.0-c)));
	// Declare B.
	std::complex<long double> B = (kappa-rho*sigma*u*i-d)/pow(sigma, 2.0)*
        (((long double)1.0 - exp(-d*tau))/((long double)1.0 - c*exp(-d*tau)));
	
	// Second compute Merton coefficients.
	// Characteristic function of normal dist.
	std::complex<long double> fj = exp(i*alpha*u - 
		pow(delta*u, 2.0)*(long double) 0.5);
	std::complex<long double> xij=exp(alpha+(long double)0.5*pow(delta, 2.0));
	
	// Declare D.
	std::complex<long double> D = ((fj-(long double)1.) - (xij-(long double)1.)*i*u)*tau;
	
	return exp(A + B*vol + D*lambda);
}

template <const int precision>
long double EuropeanOptionBates(
	long double tau, long double spot, long double vol, long double strike, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda,
    long double r, long double q) {
	// Bind the characteristic funtion.
	std::function<std::complex<long double> (std::complex<long double>)> phi = 
		std::bind(CharacteristicFunctionBates, std::placeholders::_1, 
		tau, vol, 
		kappa, theta, sigma, rho,
		alpha, delta, lambda);
	
    // Compute price.
    long double price = EuropeanOptionFourier<precision>(phi, tau, spot, strike, r, q);
    return price;
}

void CallbackEuropeanOptionBates(const alglib::real_1d_array &c, const alglib::real_1d_array &x, 
    double &func, void *ptr) {
    func = EuropeanOptionBates<5000>(
		x[0], x[1], c[0], x[2],
		c[1], c[2], c[3], c[4],
		c[5], c[6], c[7],
		x[3], x[4]);
}

boost::numeric::ublas::vector<double> CalibrateBates(boost::numeric::ublas::matrix<double> call_data, long double spot) {
	using namespace alglib;
	// The expected form of call_data is five columns with maturity, rate, dividend, strike, impvol.
    int numRows = call_data.size1();
    
	// Create arrays.
    real_2d_array x;
    x.setlength(numRows, 7);
    
    real_1d_array y;
    y.setlength(numRows);
    
    // Define parameters and bounds.
    //                      vol, kappa, theta, sigma, rho, alpha, delta, lambda.
    real_1d_array c    = "[0.0672, 2.00, 0.10,  0.2, -0.5, -0.3,  0.2,  0.2]"; // initial guess
    real_1d_array bndl = "[0.00, 0.00, 0.00,  0.0, -1.0, -inf,  0.0,  0.0]";
    real_1d_array bndu = "[+inf, +inf, +inf, +inf,  1.0, +inf, +inf, +inf]";
		
    // Fill arrays.
    long double tau, strike, r, q, impvol;
	long double deflator, F, xx, price;
	for(int i = 0; i < numRows; i++) { 
        // Extract values.
		tau = call_data(i, 0); 	  // maturity
		r = call_data(i, 1); 	  // rate
		q = call_data(i, 2); 	  // dividend
                strike = call_data(i, 3);  // strike
		impvol = call_data(i, 4);  // impvol
		
		price = EuropeanOptionBlackScholes(tau, spot, strike, impvol, r, q, 1, 1);
		
		// fill x
		x(i, 0) = tau;
		x(i, 1) = spot;
		x(i, 2) = strike;
		x(i, 3) = r;
		x(i, 4) = q;
			
		// Fill y.
        y(i) = price;
    }
    
	// Basic setup
    double epsx = 0.000001;
    ae_int_t maxits = 0;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;
    double diffstep = 0.0001;
    
    // Fit parameters.
    lsfitcreatef(x, y, c, diffstep, state);
    lsfitsetbc(state, bndl, bndu);
    lsfitsetcond(state, epsx, maxits);
    alglib::lsfitfit(state, CallbackEuropeanOptionBates);
    lsfitresults(state, info, c, rep);
    
	// Summary of results.
	// -8 : optimizer   detected  NAN/INF  in  the  target function and/or gradient
	// -7 : gradient verification failed. See LSFitSetGradientCheck() for more information.
	// -3 : inconsistent constraints
	//  2 : relative step is no more than EpsX.
	//  5 : MaxIts steps was taken.
	//  7 : stopping conditions are too stringent, further improvement is impossible.
    printf("Completion code: %d\n", int(info));
    printf("Solution: %s\n", c.tostring(2).c_str());
	
	// Calculate statistics.
	long double squared_error = 0;
	long double squared_error_short_term = 0;
	long double approxPrice, approxImpvol;
	for(int i = 0; i < numRows; i++) {
		approxPrice = EuropeanOptionBates<5000>(
			x(i, 0), x(i, 1), c[0], x(i, 2),
			c[1], c[2], c[3], c[4],
			c[5], c[6], c[7],
			x(i, 3), x(i, 4));
		approxImpvol = ImpliedVolatilityJackel(approxPrice, x(i, 0), x(i, 1), x(i, 2), x(i, 3), x(i, 4));	
		squared_error += pow(approxImpvol - call_data(i, 4), 2.);
		
		// calculate short term error
		if(x(i, 0) < 0.1) {
			squared_error_short_term += pow(approxImpvol - call_data(i, 4), 2.);
		}
	}
	squared_error = sqrt(squared_error);
	squared_error_short_term = sqrt(squared_error_short_term);
	
	printf("Euclidian norm overall: %.15Lg\n", (double) squared_error);
	printf("Euclidian norm short term: %.15Lg\n", (double) squared_error_short_term);

	// Fill return value with parameters.
	// kappa, sigma, rho, alpha, delta, lambda.
    boost::numeric::ublas::vector<double> fitted_params(9);
    fitted_params(0) = spot;   // spot
	fitted_params(1) = c[0];   // vol
	fitted_params(2) = c[1];   // kappa
	fitted_params(3) = c[2];   // theta
	fitted_params(4) = c[3];   // sigma
	fitted_params(5) = c[4];   // rho
	fitted_params(6) = c[5];   // alpha
	fitted_params(7) = c[6];   // delta
	fitted_params(8) = c[7];   // lambda
	
	return fitted_params;
}

// Template definitions
template long double EuropeanOptionBates<5000>(
	long double tau, long double spot, long double vol, long double strike, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda,
    long double r, long double q);

template long double EuropeanOptionBates<50000>(
	long double tau, long double spot, long double vol, long double strike, 
	long double kappa, long double theta, long double sigma, long double rho,
	long double alpha, long double delta, long double lambda,
    long double r, long double q);
