#include "hdp_dynamic.h"
#include <assert.h>

#define VERBOSE true

#define PERMUTE true
#define PERMUTE_LAG 10
#define TABLE_SAMPLING true
#define INF -1e50

hdp_dynamic::hdp_dynamic() {
	m_hdp_param = NULL;
	m_state = NULL;
}

hdp_dynamic::~hdp_dynamic() {
	if (m_hdp_param != NULL) {
		delete m_hdp_param;
	}
	m_hdp_param = NULL;

	if (m_state != NULL) {
		delete m_state;
	}
	m_state = NULL;
}

void hdp_dynamic::run() {
	m_state->iterate_gibbs_state(false, m_hdp_param, TABLE_SAMPLING); //init the state

	printf("starting with %d topics \n", m_state->m_num_topics);

	/// time part ...
	time_t start, current;
	time(&start);
	double dif;
	///
	for (int iter = 0; iter < m_hdp_param->m_max_iter; iter++) {
		printf("iter = %05d, ", iter);

		m_state->iterate_gibbs_state(true, m_hdp_param, TABLE_SAMPLING);

		time(&current);
		dif = difftime(current, start);

		printf("#topics = %04d, time = %8.2f\n", m_state->m_num_topics, dif);

	}
}

void hdp_dynamic::run_online() {
	m_state->iterate_gibbs_state_online(false, m_hdp_param, TABLE_SAMPLING); //init the state

	for (int iter = 0; iter < m_hdp_param->m_max_iter; iter++) {
		m_state->iterate_gibbs_state_online(true, m_hdp_param, TABLE_SAMPLING);
	}
}

void hdp_dynamic::setup_state(const corpus *c, double _eta, double _gamma,
		double _alpha, double _delta, const hdp_hyperparameter _hdp_param,
		int vocab_size) {
	free_hdp();

	m_hdp_param = new hdp_hyperparameter();
	m_hdp_param->copy_parameters(_hdp_param);
	m_state = new hdp_state_dynamic();

	m_state->set_vocab_size(vocab_size);
	m_state->setup_state_from_corpus(c);
	m_state->allocate_initial_space();
	m_state->m_eta = _eta;
	m_state->m_delta = _delta;
	/// set the pointer for the hyper parameters
	assert(m_hdp_param != NULL);
	/// use the means of gamma distribution
	if (_gamma <= 0) {
		m_state->m_gamma = m_hdp_param->m_gamma_a * m_hdp_param->m_gamma_b;
	} else {
		m_state->m_gamma = _gamma;
	}
	if (_alpha <= 0) {
		m_state->m_alpha = m_hdp_param->m_alpha_a * m_hdp_param->m_alpha_b;
	} else {
		m_state->m_alpha = _alpha;
	}
}

void hdp_dynamic::setup_state(const hdp_hyperparameter _hdp_param) {
	if (m_hdp_param != NULL) {
		delete m_hdp_param;
		m_hdp_param = NULL;
	}
	m_hdp_param = new hdp_hyperparameter();
	m_hdp_param->copy_parameters(_hdp_param);
}

void hdp_dynamic::setup_doc_info_update(const document * doc) {
	m_state->setup_doc_info_from_document(doc);
	m_state->allocate_space_for_online_updates();
}

int hdp_dynamic::get_size_vocab() const {
	if (m_state != NULL) {
		return m_state->m_size_vocab;
	} else {
		return -1;
	}
}

double hdp_dynamic::get_eta() const {
	if (m_state != NULL) {
		return m_state->m_eta;
	} else {
		return -1;
	}
}

int hdp_dynamic::get_topics_number() const {
	if (m_state != NULL) {
		return m_state->m_num_topics;
	} else {
		return -1;
	}
}

int hdp_dynamic::get_word_topic_counts(int k) const {
	if (m_state != NULL) {
		return m_state->m_word_counts_by_z[k];
	} else {
		return -1;
	}
}

int hdp_dynamic::get_word_topic_counts(int k, int w) const {
	if (m_state != NULL) {
		return m_state->m_word_counts_by_zw[k][w];
	} else {
		return -1;
	}
}

void hdp_dynamic::load(char * model_path, int sampler_num) {
	free_hdp();
	m_state = new hdp_state_dynamic();
	m_state->load_state_ex(model_path, sampler_num);
}

void hdp_dynamic::save(char * model_path, int sample_num) {
	m_state->save_state_ex(model_path, sample_num);
}

void hdp_dynamic::free_doc_info() {
	m_state->free_doc_info();
}

void hdp_dynamic::update_state_after_online_update() {
	m_state->update_state_after_online_update();
}

void hdp_dynamic::free_hdp() {
	if (m_hdp_param != NULL) {
		delete m_hdp_param;
	}
	if (m_state != NULL) {
		delete m_state;
	}
}

void hdp_dynamic::save_last_document_info_for_online_update() {
	m_state->save_last_document_info_for_online_update();
}
