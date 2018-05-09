#include "hdp.h"
#include "likelihood_computer.h"
#include "utils.h"
#include "corpus.h"
#include "state.h"
#include <map>
#include <iostream>
#include <vector>
using namespace std;

double likelihood_computer::compute_harmonic_mean_predictive_likelihood(
		const hdp* const * posterior_samplers,
		const counts* const * before_update_counts, int m_posterior_samplers) {
	double likelihood = 0;
	if (posterior_samplers == NULL) {
		return likelihood;
	}
	int m_size_vocab = posterior_samplers[0]->get_size_vocab();
	double m_eta = posterior_samplers[0]->get_eta();

	vector<double> log_likelihood;
	log_likelihood.assign(m_posterior_samplers, 0);
	for (int posterior_sampler_id = 0;
			posterior_sampler_id < m_posterior_samplers;
			++posterior_sampler_id) {

		int m_old_topics =
				before_update_counts[posterior_sampler_id]->m_word_counts_by_z.size();
		int m_new_topics =
				posterior_samplers[posterior_sampler_id]->get_topics_number();

		for (int k = 0; k < m_old_topics; ++k) {
			log_likelihood[posterior_sampler_id] +=
					lgamma(
							m_size_vocab * m_eta
									+ before_update_counts[posterior_sampler_id]->m_word_counts_by_z[k])
							- lgamma(
									m_size_vocab * m_eta
											+ posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
													k));
			for (int w = 0; w < m_size_vocab; ++w) {
				int cur_word_topic_count =
						posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
								k, w);
				int cur_word_topic_count_before_update =
						before_update_counts[posterior_sampler_id]->m_word_counts_by_zw[k][w];
				if (cur_word_topic_count
						!= cur_word_topic_count_before_update) {
					log_likelihood[posterior_sampler_id] += lgamma(
							cur_word_topic_count + m_eta)
							- lgamma(
									cur_word_topic_count_before_update + m_eta);
				}
			}
		}
		for (int k = m_old_topics; k < m_new_topics; ++k) {
			log_likelihood[posterior_sampler_id] +=
					lgamma(m_size_vocab * m_eta)
							- lgamma(
									m_size_vocab * m_eta
											+ posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
													k));
			for (int w = 0; w < m_size_vocab; ++w) {
				if (posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
						k, w) > 0) {
					log_likelihood[posterior_sampler_id] +=
							lgamma(
									posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
											k, w) + m_eta) - lgamma(m_eta);
				}
			}
		}
	}

	invert_element_sign(log_likelihood);

	likelihood = -log_sum_exp_trick(log_likelihood) + log(m_posterior_samplers);

	return likelihood;
}

double likelihood_computer::compute_harmonic_mean_predictive_likelihood(
		const hdp_dynamic* const * posterior_samplers,
		const counts* const * before_update_counts, int m_posterior_samplers) {
	double likelihood = 0;
	if (posterior_samplers == NULL) {
		return likelihood;
	}
	int m_size_vocab = posterior_samplers[0]->get_size_vocab();
	double m_eta = posterior_samplers[0]->get_eta();

	vector<double> log_likelihood;
	log_likelihood.assign(m_posterior_samplers, 0);
	for (int posterior_sampler_id = 0;
			posterior_sampler_id < m_posterior_samplers;
			++posterior_sampler_id) {

		int m_old_topics =
				before_update_counts[posterior_sampler_id]->m_word_counts_by_z.size();
		int m_new_topics =
				posterior_samplers[posterior_sampler_id]->get_topics_number();

		for (int k = 0; k < m_old_topics; ++k) {
			log_likelihood[posterior_sampler_id] +=
					lgamma(
							m_size_vocab * m_eta
									+ before_update_counts[posterior_sampler_id]->m_word_counts_by_z[k])
							- lgamma(
									m_size_vocab * m_eta
											+ posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
													k));
			for (int w = 0; w < m_size_vocab; ++w) {
				int cur_word_topic_count =
						posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
								k, w);
				int cur_word_topic_count_before_update =
						before_update_counts[posterior_sampler_id]->m_word_counts_by_zw[k][w];
				if (cur_word_topic_count
						!= cur_word_topic_count_before_update) {
					log_likelihood[posterior_sampler_id] += lgamma(
							cur_word_topic_count + m_eta)
							- lgamma(
									cur_word_topic_count_before_update + m_eta);
				}
			}
		}
		for (int k = m_old_topics; k < m_new_topics; ++k) {
			log_likelihood[posterior_sampler_id] +=
					lgamma(m_size_vocab * m_eta)
							- lgamma(
									m_size_vocab * m_eta
											+ posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
													k));
			for (int w = 0; w < m_size_vocab; ++w) {
				if (posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
						k, w) > 0) {
					log_likelihood[posterior_sampler_id] +=
							lgamma(
									posterior_samplers[posterior_sampler_id]->get_word_topic_counts(
											k, w) + m_eta) - lgamma(m_eta);
				}
			}
		}
	}

	invert_element_sign(log_likelihood);

	likelihood = -log_sum_exp_trick(log_likelihood) + log(m_posterior_samplers);

	return likelihood;
}
