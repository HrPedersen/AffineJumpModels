#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <functional>

// Functions
long double dnorm(long double x) {
    return (long double)1.0/((long double)2.0*M_PI)*exp(-(long double)0.5*pow(x, 2.0));
}

long double pnorm(long double x) {
   return (long double)0.5 * erfc(-x * M_SQRT1_2);
}

long double EuropeanOptionBlackScholes(long double tau, long double spot, 
	long double strike, long double sigma, 
    double r, double q, int greek, int option) {
    // Compute dummies.
    long double d1 = (log(spot/strike)+ ((r-q)+0.5*pow(sigma, 2.0))*tau)/
        (sigma*sqrt(tau));
    long double d2 = d1 - sigma*sqrt(tau);
    
    // Calculate price
    long double price;
    if(greek == 1 & option == 1) {
		price = spot*exp(-q*tau)*pnorm(d1)-strike*exp(-r*tau)*pnorm(d2);
		return price;
	} else if(greek == 1 & option == 2) {
		price = spot*exp(-q*tau)*pnorm(d1)-strike*exp(-r*tau)*pnorm(d2)-
			spot*exp(-q*tau)+strike*exp(-r*tau);
		return price;
	}
	
	
    // Calculate first order greeks.
    long double delta = exp(-q*tau)*pnorm(d1);
    long double vega = exp(-q*tau)*spot*dnorm(d1)*sqrt(tau);
    
    // Calculate second order greeks.
    long double gamma = exp(-q*tau)*dnorm(d1)/(spot*sigma*sqrt(tau));
    long double vanna = vega/spot*(1-d1/(sigma*sqrt(tau)));
    long double volga = vega*d1*d2/sigma;
    
    // Return proper greek.
    if(greek == 2) {
        return delta;
    } else if(greek == 3) {
        return vega;
    } else if(greek == 4) {
        return gamma;
    } else if(greek == 5) {
        return vanna;
    } else if(greek == 6) {
        return volga;
    }
}

long double cdf(long double x) {
	return 0.5 * erfc(-x * M_SQRT1_2);
}

long double inv_cdf(long double p) {
    // https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/
    // Coefficients in rational approximations.
    long double q, x, r; // output.
	
	long double a1 = -3.969683028665376e+01;
    long double a2 =  2.209460984245205e+02;
    long double a3 = -2.759285104469687e+02;
    long double a4 =  1.383577518672690e+02;
    long double a5 = -3.066479806614716e+01;
    long double a6 =  2.506628277459239e+00;
	
    long double b1 = -5.447609879822406e+01;
    long double b2 =  1.615858368580409e+02;
    long double b3 = -1.556989798598866e+02;
    long double b4 =  6.680131188771972e+01;
    long double b5 = -1.328068155288572e+01;
	
    long double c1 = -7.784894002430293e-03;
    long double c2 = -3.223964580411365e-01;
    long double c3 = -2.400758277161838e+00;
    long double c4 = -2.549732539343734e+00;
    long double c5 =  4.374664141464968e+00;
    long double c6 =  2.938163982698783e+00;
	
    long double d1 =  7.784695709041462e-03;
    long double d2 =  3.224671290700398e-01;
    long double d3 =  2.445134137142996e+00;
    long double d4 =  3.754408661907416e+00;

    // Define break-points.
	long double p_low = 0.02425;
	long double p_high = 1.0 - p_low;

	// Rational approximation for lower region.
    if (0.0 < p < p_low) {
		q = sqrt(-2.0*log(p));
		x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
	}

    // Rational approximation for central region.
	if (p_low <= p <= p_high) {
		q = p - 0.5;
		r = q*q;
		x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
           (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0);
	}

    // Rational approximation for upper region.
	if (p_high < p < 1.0) {
		q = sqrt(-2.0*log(1.0-p));
		x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0);
	}
	
	return x;
}

long double BlackScholesNormalised(long double x, long double sigma, 
	long double deflator, long double theta) {
	throw std::runtime_error("Not implemented");
	
	// For simplicity.
	long double dp = x/sigma + sigma/(long double)2.0;
	long double dm = x/sigma - sigma/(long double)2.0;
	
	return deflator*(theta*exp(x/(long double)2.0)*cdf(theta*dp) - 
		theta*exp(-x/(long double)2.0)*cdf(theta*dm));
}
