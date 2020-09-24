#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "def.h"

namespace mips {

// -----------------------------------------------------------------------------
inline float uniform(				// r.v. from Uniform(min, max)
	float min,							// min value
	float max)							// max value
{
	int   num  = rand();
	float base = (float) RAND_MAX - 1.0F;
	float frac = ((float) num) / base;

	return (max - min) * frac + min;
}

// -----------------------------------------------------------------------------
//	Given a mean and a standard deviation, gaussian generates a normally 
//		distributed random number.
//
//	Algorithm:  Polar Method, p.104, Knuth, vol. 2
// -----------------------------------------------------------------------------
float gaussian(				// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma);						// std value

// -----------------------------------------------------------------------------
inline float normal_pdf(			// pdf of Guassian(mean, std)
	float x,							// variable
	float u,							// mean
	float sigma)						// standard error
{
	float ret = exp(-(x - u) * (x - u) / (2.0f * sigma * sigma));
	ret /= sigma * sqrt(2.0f * PI);
	return ret;
}

// -----------------------------------------------------------------------------
float normal_cdf(					// cdf of N(0, 1) in (-inf, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float new_cdf(						// cdf of N(0, 1) in [-x, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

} // end namespace mips
