#ifndef __RANDOM_H
#define __RANDOM_H

// -----------------------------------------------------------------------------
float uniform(						// r.v. from Uniform(min, max)
	float min,							// min value
	float max);							// max value

// -----------------------------------------------------------------------------
float gaussian(						// r.v. from Gaussian(mean, sigma)
	float mean,							// mean value
	float sigma);						// std value

// -----------------------------------------------------------------------------
float normal_pdf(					// pdf of Guassian(mean, std)
	float x,							// variable
	float u,							// mean
	float sigma);						// standard error

// -----------------------------------------------------------------------------
float normal_cdf(					// cdf of N(0, 1) in (-inf, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

// -----------------------------------------------------------------------------
float new_cdf(						// cdf of N(0, 1) in [-x, x]
	float x,							// integral border
	float step = 0.001f);				// step increment

#endif // __RANDOM_H
