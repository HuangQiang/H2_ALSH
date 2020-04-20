#include "random.h"

// -----------------------------------------------------------------------------
//	Given a mean and a standard deviation, gaussian generates a normally 
//		distributed random number.
//
//	Algorithm:  Polar Method, p.104, Knuth, vol. 2
// -----------------------------------------------------------------------------
float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma)						// std value
{
	float v1 = -1.0f;
    float v2 = -1.0f;
	float s  = -1.0f;
	float x  = -1.0f;

	do {
		v1 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		v2 = 2.0F * uniform(0.0F, 1.0F) - 1.0F;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.0F);
	x = v1 * sqrt(-2.0F * log (s) / s);

	x = x * sigma + mean; 			// x is distributed from N(0, 1)
	return x;
}

// -----------------------------------------------------------------------------
float normal_cdf(					// cdf of N(0, 1) in range (-inf, x]
	float _x,							// integral border
	float _step)						// step increment
{
	float ret = 0.0;
	for (float i = -10.0; i < _x; i += _step) {
		ret += _step * normal_pdf(i, 0.0f, 1.0f);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float new_cdf(						// cdf of N(0, 1) in range [-x, x]
	float x,							// integral border
	float step)							// step increment
{
	float result = 0.0f;
	for (float i = -x; i <= x; i += step) {
		result += step * normal_pdf(i, 0.0f, 1.0f);
	}
	return result;
}
