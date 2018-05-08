#ifndef LIKELIHOOD_COMPUTER_H
#define LIKELIHOOD_COMPUTER_H

#include "hdp.h"
#include "corpus.h"

class likelihood_computer {
public:
	double compute_harmonic_mean_predictive_likelihood(const hdp* const * posterior_samplers,
													   const counts* const* before_update_counts, 
													   int m_posterior_samplers);
	double compute_harmonic_mean_predictive_likelihood(const hdp_dynamic* const * posterior_samplers, 
													   const counts* const* before_update_counts, 
													   int m_posterior_samplers);
};


#endif