
#ifdef __cplusplus
extern "C" {
#endif

#ifndef GSL_GAUSS
#define GLS_GAUSS

#include <math.h>
#include <gsl_math.h>
#include <gsl_cdf.h>

double gsl_cdf_ugaussian_P (const double x);
double gsl_cdf_ugaussian_Q (const double x);
double gsl_cdf_gaussian_P (const double x, const double sigma);
double gsl_cdf_gaussian_Q (const double x, const double sigma);

#endif // GLS_GAUSS

#ifdef __cplusplus
}
#endif
