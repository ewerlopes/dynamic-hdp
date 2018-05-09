#include "state.h"
#include "utils.h"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

#define SMALL_GIBBS_MAX_ITER 20
#define VERBOSE false
#define INIT_SIZE 50
#define INIT_ONLINE_DOCS 100
#define INF -1e50
#define EPS 0.00000000001

void hdp_hyperparameter::copy_parameters(
		const hdp_hyperparameter _input_hyperparam) {
	m_gamma_a = _input_hyperparam.m_gamma_a;
	m_gamma_b = _input_hyperparam.m_gamma_b;
	m_alpha_a = _input_hyperparam.m_alpha_a;
	m_alpha_b = _input_hyperparam.m_alpha_b;
	m_max_iter = _input_hyperparam.m_max_iter;
}

doc_state::doc_state() {
	m_doc_id = -1; // document id
	m_doc_length = 0;  // document length
	m_num_tables = 0;  // number of tables in this document
	m_words = NULL;

	m_table_to_topic.clear(); // for a doc, translate its table index to topic index
	m_word_counts_by_t.clear(); // word counts for each table
}

doc_state::~doc_state() {
	free_doc_state();
}

void doc_state::setup_state_from_doc(const document * doc) {
	m_doc_id = doc->id;
	m_doc_length = doc->total;

	int word, count;
	m_words = new word_info[doc->total];
	int m = 0;
	for (int n = 0; n < doc->length; n++) {
		word = doc->words[n];
		count = doc->counts[n];
		for (int j = 0; j < count; j++) {
			m_words[m].m_word_index = word;
			m_words[m].m_table_assignment = -1;
			//m_words[m].m_topic_assignment = -1;
			m++;
		}
	}
	//allocate some space
	m_table_to_topic.resize(INIT_SIZE, -1);
	m_word_counts_by_t.resize(INIT_SIZE, 0);
}

void doc_state::free_doc_state() {
	m_table_to_topic.clear(); // for a doc, translate its table index to topic index
	m_word_counts_by_t.clear(); // word counts for each table

	delete[] m_words;
	m_words = NULL;
}

counts::counts() {
	m_num_words = 0;
	m_total_num_tables = 0;
	m_num_tables_by_z.clear();
	m_word_counts_by_z.clear();
	m_word_counts_by_zw.clear();
}

counts::counts(int _m_num_words, int _m_total_num_tables,
		const int_vec& _m_num_tables_by_z, const int_vec& _m_word_counts_by_z,
		const vector<int*>& _m_word_counts_by_zw) {
	m_num_words = _m_num_words;
	m_total_num_tables = _m_total_num_tables;
	m_num_tables_by_z = _m_num_tables_by_z;
	m_word_counts_by_z = _m_word_counts_by_z;
	int m_num_topics = m_num_tables_by_z.size();
	m_word_counts_by_zw.resize(m_num_topics, NULL);
	for (int k = 0; k < m_num_topics; ++k) {
		int* p = new int[m_num_words];
		memset(p, 0, sizeof(int) * m_num_words);
		m_word_counts_by_zw[k] = p;
		copy(_m_word_counts_by_zw[k], _m_word_counts_by_zw[k] + m_num_words,
				m_word_counts_by_zw[k]);
//        for (int i = 0; i < m_num_words; ++i) {
//            m_word_counts_by_zw[k][i] = _m_word_counts_by_zw[k][i];
//        }
	}
}

counts::~counts() {
	m_num_words = 0;
	m_total_num_tables = 0;
	m_word_counts_by_z.clear();
	m_num_tables_by_z.clear();

	for (int i = 0; i < m_word_counts_by_zw.size(); ++i) {
		if (!(m_word_counts_by_zw[i] == NULL)) {
			delete[] m_word_counts_by_zw[i];
		}
	}
	m_word_counts_by_zw.clear();
}

void counts::set_counts(int _m_num_words, int _m_total_num_tables,
		const int_vec& _m_num_tables_by_z, const int_vec& _m_word_counts_by_z,
		const vector<int*>& _m_word_counts_by_zw) {
	if ((m_num_words == _m_num_words)
			&& (m_num_tables_by_z.size() == _m_num_tables_by_z.size())) {
		// just change the values
		m_total_num_tables = _m_total_num_tables;
		m_num_tables_by_z = _m_num_tables_by_z;
		m_word_counts_by_z = _m_word_counts_by_z;
		int m_num_topics = m_num_tables_by_z.size();
		for (int k = 0; k < m_num_topics; ++k) {
			for (int i = 0; i < m_num_words; ++i) {
				m_word_counts_by_zw[k][i] = _m_word_counts_by_zw[k][i];
			}
		}
	} else {
		// free and allocate memory as well as set new values
		m_num_tables_by_z.clear();
		m_word_counts_by_z.clear();

		int size = m_word_counts_by_zw.size();
		for (int k = 0; k < size; ++k) {
			if (!(m_word_counts_by_zw[k] == NULL))
				delete[] m_word_counts_by_zw[k];
		}
		m_word_counts_by_zw.clear();

		m_num_words = _m_num_words;
		m_total_num_tables = _m_total_num_tables;
		m_num_tables_by_z = _m_num_tables_by_z;
		m_word_counts_by_z = _m_word_counts_by_z;
		int m_num_topics = m_num_tables_by_z.size();
		m_word_counts_by_zw.resize(m_num_topics, NULL);
		for (int k = 0; k < m_num_topics; ++k) {
			int* p = new int[m_num_words];
			memset(p, 0, sizeof(int) * m_num_words);
			m_word_counts_by_zw[k] = p;
			copy(_m_word_counts_by_zw[k], _m_word_counts_by_zw[k] + m_num_words,
					m_word_counts_by_zw[k]);
//			for (int i = 0; i < m_num_words; ++i) {
//				m_word_counts_by_zw[k][i] = _m_word_counts_by_zw[k][i];
//			}
		}
	}
}

word_counts::word_counts() {
	m_word_counts_by_z.clear();
	m_word_counts_by_zw.clear();
}

word_counts::~word_counts() {
	m_word_counts_by_z.clear();
	for (int k = 0; k < m_word_counts_by_zw.size(); ++k) {
		if (m_word_counts_by_zw[k] != NULL) {
			delete[] m_word_counts_by_zw[k];
		}
	}
	m_word_counts_by_zw.clear();
}

hdp_state::hdp_state() {
	m_doc_states = NULL;
	m_size_vocab = 0;
	m_total_words = 0;
	;
	m_num_docs = 0;
	m_num_topics = 0;
	m_total_num_tables = 0;
	m_num_tables_by_z.clear();
	m_word_counts_by_z.clear();
	m_word_counts_by_zd.clear();
	m_word_counts_by_zw.clear();

	m_word_counts_by_zd_online.clear();
}

hdp_state::~hdp_state() {
	free_state();
}

void hdp_state::set_vocab_size(int _vocab_size) {
	m_size_vocab = _vocab_size;
}

void hdp_state::setup_state_from_corpus(const corpus * c) {
	m_size_vocab = max(c->size_vocab, m_size_vocab);
	m_total_words += c->total_words;
	m_num_docs = c->num_docs;
	m_doc_states = new doc_state *[m_num_docs];

	for (unsigned int d = 0; d < c->docs.size(); d++) {
		document * doc = c->docs[d];
		doc_state * d_state = new doc_state();
		m_doc_states[d] = d_state;
		d_state->setup_state_from_doc(doc);
	}
}

void hdp_state::setup_doc_info_from_document(const document * doc) {
	int proposed_vocab_size = doc->get_max_word_id() + 1;
	m_size_vocab = max(m_size_vocab, proposed_vocab_size);
	m_total_words += doc->total;

	if (!(m_doc_states == NULL)) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
			d_state = NULL;
		}
		delete[] m_doc_states;
		m_doc_states = NULL;
	}
	m_num_docs = 1;
	m_doc_states = new doc_state *[1];

	doc_state * d_state = new doc_state();
	m_doc_states[0] = d_state;
	d_state->setup_state_from_doc(doc);
}

void hdp_state::allocate_space_for_online_updates() {
	m_word_counts_by_zd_online.resize(m_num_topics, 0);

	m_num_topics_before_update = m_num_topics;
}

void hdp_state::allocate_initial_space() {
	// training
	if (m_num_tables_by_z.size() == 0) {
		m_num_tables_by_z.resize(INIT_SIZE, 0);
		m_word_counts_by_z.resize(INIT_SIZE, 0);
		m_word_counts_by_zd.resize(INIT_SIZE, NULL);
		m_word_counts_by_zw.resize(INIT_SIZE, NULL);
		int * p = NULL;
		for (int k = 0; k < INIT_SIZE; k++) {
			p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_word_counts_by_zd[k] = p;

			p = new int[m_size_vocab];
			memset(p, 0, sizeof(int) * m_size_vocab);
			m_word_counts_by_zw[k] = p;
		}
	} else // testing
	{
		while ((int) m_num_tables_by_z.size() < m_num_topics + 1) {
			m_num_tables_by_z.push_back(0);
			m_word_counts_by_z.push_back(0);

			int* p = new int[m_size_vocab];
			memset(p, 0, sizeof(int) * m_size_vocab);
			m_word_counts_by_zw.push_back(p);
		}
		while ((int) m_word_counts_by_zd.size() < m_num_topics + 1) {
			int* p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_word_counts_by_zd.push_back(p);
		}
	}
}

void hdp_state::free_state() {
	if (m_doc_states != NULL) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
		}
		delete[] m_doc_states;
	}
	m_doc_states = NULL;

	m_size_vocab = 0;
	m_total_words = 0;
	;
	m_num_docs = 0;
	m_num_topics = 0;
	m_num_topics_before_update = 0;
	m_total_num_tables = 0;

	m_num_tables_by_z.clear();
	m_word_counts_by_z.clear();

	free_vec_ptr(m_word_counts_by_zd);
	free_vec_ptr(m_word_counts_by_zw);

	m_word_counts_by_zd_online.clear();
}

void hdp_state::free_doc_info() {
	if (m_doc_states != NULL) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
		}
		delete[] m_doc_states;
	}
	m_doc_states = NULL;

	m_num_docs = 0;

	for (int k = 0; k < m_word_counts_by_zd.size(); ++k) {
		delete[] m_word_counts_by_zd[k];
	}
	m_word_counts_by_zd.clear();
}

void hdp_state::free_counts_update() {
	m_word_counts_by_zd_online.clear();
	m_num_topics_before_update = 0;
}

void hdp_state::free_online_info() {
	free_doc_info();
	free_counts_update();
}

void hdp_state::iterate_gibbs_state(bool remove, bool permute,
		hdp_hyperparameter * hdp_hyperparam, bool table_sampling) {
	if (permute) // permuate data
	{
		rshuffle(m_doc_states, m_num_docs, sizeof(doc_state*));
		for (int j = 0; j < m_num_docs; j++)
			rshuffle(m_doc_states[j]->m_words, m_doc_states[j]->m_doc_length,
					sizeof(word_info));
	}

	double_vec q;
	double_vec f;
	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		d_state = m_doc_states[j];
		for (int i = 0; i < d_state->m_doc_length; i++) {
			sample_word_assignment(d_state, i, remove, q, f);
		}
		if (table_sampling)
			sample_tables(d_state, q, f);

	}
	compact_hdp_state();

}

void hdp_state::iterate_gibbs_state_online(bool remove,
		hdp_hyperparameter* hdp_hyperparam, bool table_sampling) {
	double_vec q;
	double_vec f;
	doc_state* d_state = m_doc_states[0];

	for (int i = 0; i < d_state->m_doc_length; i++) {
		sample_word_assignment_online(d_state, i, remove, q, f);
	}
	if (table_sampling)
		sample_tables_online(d_state, q, f);

	compact_hdp_state_online();
}

void hdp_state::compact_doc_state_online(doc_state* d_state, int* k_to_new_k) {
	int num_tables_old = d_state->m_num_tables;
	int* t_to_new_t = new int[num_tables_old];

	int t, new_t, k, w;
	for (t = 0, new_t = 0; t < num_tables_old; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			t_to_new_t[t] = new_t;
			k = d_state->m_table_to_topic[t];
			if (k >= m_num_topics_before_update) {
				d_state->m_table_to_topic[new_t] = k_to_new_k[k
						- m_num_topics_before_update];
			} else {
				d_state->m_table_to_topic[new_t] = k;
			}
			swap_vec_element(d_state->m_word_counts_by_t, new_t, t);
			new_t++;
		} else
			d_state->m_table_to_topic[t] = -1;
	}
	d_state->m_num_tables = new_t;

	for (int i = 0; i < d_state->m_doc_length; i++) {
		t = d_state->m_words[i].m_table_assignment;
		new_t = t_to_new_t[t];
		d_state->m_words[i].m_table_assignment = new_t;
		w = d_state->m_words[i].m_word_index;
	}

	delete[] t_to_new_t;
}

void hdp_state::compact_doc_state(doc_state* d_state, int* k_to_new_k) {
	int num_tables_old = d_state->m_num_tables;
	int* t_to_new_t = new int[num_tables_old];

	int t, new_t, k, w;
	for (t = 0, new_t = 0; t < num_tables_old; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			t_to_new_t[t] = new_t;
			k = d_state->m_table_to_topic[t];
			d_state->m_table_to_topic[new_t] = k_to_new_k[k];
			swap_vec_element(d_state->m_word_counts_by_t, new_t, t);
			new_t++;
		} else
			d_state->m_table_to_topic[t] = -1;
	}
	d_state->m_num_tables = new_t;

	for (int i = 0; i < d_state->m_doc_length; i++) {
		t = d_state->m_words[i].m_table_assignment;
		new_t = t_to_new_t[t];
		d_state->m_words[i].m_table_assignment = new_t;
		w = d_state->m_words[i].m_word_index;
	}

	delete[] t_to_new_t;
}

//compress the unused tables and components
void hdp_state::compact_hdp_state() {
	int num_topics_old = m_num_topics;
	int* k_to_new_k = new int[num_topics_old];
	int k, new_k;
	for (k = 0, new_k = 0; k < num_topics_old; k++) {
		if (m_word_counts_by_z[k] > 0) {
			k_to_new_k[k] = new_k;
			swap_vec_element(m_word_counts_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_z, new_k, k);
			swap_vec_element(m_word_counts_by_zd, new_k, k);
			swap_vec_element(m_word_counts_by_zw, new_k, k);
			new_k++;
		}
	}
	m_num_topics = new_k;

	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		d_state = m_doc_states[j];
		compact_doc_state(d_state, k_to_new_k);
	}

	delete[] k_to_new_k;
}

//compress the unused tables and components
void hdp_state::compact_hdp_state_online() {
	int num_topics_old = m_num_topics;
	int* k_to_new_k = new int[m_num_topics - m_num_topics_before_update];
	int k, new_k;
	for (k = m_num_topics_before_update, new_k = m_num_topics_before_update;
			k < num_topics_old; k++) {
		if (m_word_counts_by_z[k] > 0) {
			k_to_new_k[k - m_num_topics_before_update] = new_k;
			swap_vec_element(m_word_counts_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_z, new_k, k);
			swap_vec_element(m_word_counts_by_zd_online, new_k, k);
			swap_vec_element(m_word_counts_by_zw, new_k, k);
			new_k++;
		}
	}
	m_num_topics = new_k;

	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		d_state = m_doc_states[j];
		compact_doc_state_online(d_state, k_to_new_k);
	}

	delete[] k_to_new_k;
}

void hdp_state::sample_tables_online(doc_state* d_state, double_vec & q,
		double_vec & f) {
	vector<int*> words_by_t;
	words_by_t.resize(d_state->m_num_tables, NULL);
	int* p = NULL;
	int t, word;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			p = new int[d_state->m_word_counts_by_t[t]];
			words_by_t[t] = p;
		}
	}

	int* i_by_t = new int[d_state->m_num_tables];
	memset(i_by_t, 0, sizeof(int) * d_state->m_num_tables);

	for (int i = 0; i < d_state->m_doc_length; i++) {
		word = d_state->m_words[i].m_word_index;
		t = d_state->m_words[i].m_table_assignment;
		words_by_t[t][i_by_t[t]] = i; // save the id for later use
		i_by_t[t]++;
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0)
			sample_table_assignment_online(d_state, t, words_by_t[t], q, f);
		//eles no needs, since there is no data there
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		p = words_by_t[t];
		delete[] p;
	}
	delete[] i_by_t;
}

void hdp_state::sample_tables(doc_state* d_state, double_vec & q,
		double_vec & f) {
	vector<int*> words_by_t;
	words_by_t.resize(d_state->m_num_tables, NULL);
	int* p = NULL;
	int t, word;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			p = new int[d_state->m_word_counts_by_t[t]];
			words_by_t[t] = p;
		}
	}

	int* i_by_t = new int[d_state->m_num_tables];
	memset(i_by_t, 0, sizeof(int) * d_state->m_num_tables);

	for (int i = 0; i < d_state->m_doc_length; i++) {
		word = d_state->m_words[i].m_word_index;
		t = d_state->m_words[i].m_table_assignment;
		words_by_t[t][i_by_t[t]] = i; // save the id for later use
		i_by_t[t]++;
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0)
			sample_table_assignment(d_state, t, words_by_t[t], q, f);
		//eles no needs, since there is no data there
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		p = words_by_t[t];
		delete[] p;
	}
	delete[] i_by_t;
}

void hdp_state::sample_table_assignment_online(doc_state* d_state, int t,
		int* words, double_vec & q, double_vec & f) {
	//number of tables won't change at all
	int i, w, k, m, k_old, d;
	int* counts = new int[m_size_vocab];
	int* counts_copy = new int[m_size_vocab];
	memset(counts_copy, 0, sizeof(int) * m_size_vocab);
	int count_sum = d_state->m_word_counts_by_t[t];

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		counts_copy[w]++;
	}
	memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);

	// compute the the log prob of being at a new cluster
	double f_new = lgamma(m_size_vocab * m_eta)
			- lgamma(count_sum + m_size_vocab * m_eta);

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		if (counts[w] > 0) {
			f_new += lgamma(counts[w] + m_eta) - lgamma(m_eta);
			counts[w] = 0;
		}
	}

	if ((int) q.size() < m_num_topics + 1)
		q.resize(2 * m_num_topics + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	q[m_num_topics] = log(m_gamma) + f_new;

	k_old = d_state->m_table_to_topic[t];

	for (k = 0; k < m_num_topics; k++) {
		if (k == k_old) {
			f[k] = lgamma(
					m_size_vocab * m_eta + m_word_counts_by_z[k] - count_sum)
					- lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k]);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(m_eta + m_word_counts_by_zw[k][w])
							- lgamma(
									m_eta + m_word_counts_by_zw[k][w]
											- counts[w]);
					counts[w] = 0;
				}
			}
			if (m_num_tables_by_z[k] == 1)
				q[k] = INF; // make it extremely small as log(0)
			else
				q[k] = log(m_num_tables_by_z[k] - 1) + f[k];
		} else {
			f[k] = lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k])
					- lgamma(
							m_size_vocab * m_eta + m_word_counts_by_z[k]
									+ count_sum);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(
							m_eta + m_word_counts_by_zw[k][w] + counts[w])
							- lgamma(m_eta + m_word_counts_by_zw[k][w]);
					counts[w] = 0;
				}
			}
			q[k] = log(m_num_tables_by_z[k]) + f[k];
		}
	}
	//normalizing in log space for sampling
	log_normalize(q, m_num_topics + 1);
	q[0] = exp(q[0]);
	double total_q = q[0];
	for (k = 1; k < m_num_topics + 1; k++) {
		total_q += exp(q[k]);
		q[k] = total_q;
	}

	double u = runiform() * total_q;
	for (k = 0; k < m_num_topics + 1; k++)
		if (u < q[k])
			break;

	if (k != k_old) // status doesn't change, but k could change
			{
		d = d_state->m_doc_id;

		/// reassign the topic to current table
		d_state->m_table_to_topic[t] = k;

		/// update the statistics by removing the table t from topic k_old
		m_num_tables_by_z[k_old]--;
		m_word_counts_by_z[k_old] -= count_sum;
		m_word_counts_by_zd_online[k_old] -= count_sum;

		if (k == m_num_topics) // a new topic is created
				{
			m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_z.size() < m_num_topics + 1) {
				m_num_tables_by_z.push_back(0);
				m_word_counts_by_z.push_back(0);

				int * p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
			if ((int) m_word_counts_by_zd_online.size() < m_num_topics) {
				m_word_counts_by_zd_online.push_back(0);
			}
		}

		/// update the statistics by adding the table t to topic k
		m_num_tables_by_z[k]++;
		m_word_counts_by_z[k] += count_sum;
		m_word_counts_by_zd_online[k] += count_sum;

		for (int m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
			i = words[m];
			w = d_state->m_words[i].m_word_index;
			m_word_counts_by_zw[k_old][w]--;
			m_word_counts_by_zw[k][w]++;
		}
	}
	delete[] counts;
	delete[] counts_copy;
}

void hdp_state::sample_table_assignment(doc_state* d_state, int t, int * words,
		double_vec & q, double_vec & f) {
	//number of tables won't change at all
	int i, w, k, m, k_old, d;
	int* counts = new int[m_size_vocab];
	int* counts_copy = new int[m_size_vocab];
	memset(counts_copy, 0, sizeof(int) * m_size_vocab);
	int count_sum = d_state->m_word_counts_by_t[t];

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		counts_copy[w]++;
	}
	memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);

	// compute the the log prob of being at a new cluster
	double f_new = lgamma(m_size_vocab * m_eta)
			- lgamma(count_sum + m_size_vocab * m_eta);

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		if (counts[w] > 0) {
			f_new += lgamma(counts[w] + m_eta) - lgamma(m_eta);
			counts[w] = 0;
		}
	}

	if ((int) q.size() < m_num_topics + 1)
		q.resize(2 * m_num_topics + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	q[m_num_topics] = log(m_gamma) + f_new;

	k_old = d_state->m_table_to_topic[t];

	for (k = 0; k < m_num_topics; k++) {
		if (k == k_old) {
			f[k] = lgamma(
					m_size_vocab * m_eta + m_word_counts_by_z[k] - count_sum)
					- lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k]);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(m_eta + m_word_counts_by_zw[k][w])
							- lgamma(
									m_eta + m_word_counts_by_zw[k][w]
											- counts[w]);
					counts[w] = 0;
				}
			}
			if (m_num_tables_by_z[k] == 1)
				q[k] = INF; // make it extremely small as log(0)
			else
				q[k] = log(m_num_tables_by_z[k] - 1) + f[k];
		} else {
			f[k] = lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k])
					- lgamma(
							m_size_vocab * m_eta + m_word_counts_by_z[k]
									+ count_sum);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(
							m_eta + m_word_counts_by_zw[k][w] + counts[w])
							- lgamma(m_eta + m_word_counts_by_zw[k][w]);
					counts[w] = 0;
				}
			}
			q[k] = log(m_num_tables_by_z[k]) + f[k];
		}
	}
	//normalizing in log space for sampling
	log_normalize(q, m_num_topics + 1);
	q[0] = exp(q[0]);
	double total_q = q[0];
	for (k = 1; k < m_num_topics + 1; k++) {
		total_q += exp(q[k]);
		q[k] = total_q;
	}

	double u = runiform() * total_q;
	for (k = 0; k < m_num_topics + 1; k++)
		if (u < q[k])
			break;

	if (k != k_old) // status doesn't change, but k could change
			{
		d = d_state->m_doc_id;

		/// reassign the topic to current table
		d_state->m_table_to_topic[t] = k;

		/// update the statistics by removing the table t from topic k_old
		m_num_tables_by_z[k_old]--;
		m_word_counts_by_z[k_old] -= count_sum;
		m_word_counts_by_zd[k_old][d] -= count_sum;

		/// update the statistics by adding the table t to topic k
		m_num_tables_by_z[k]++;
		m_word_counts_by_z[k] += count_sum;
		m_word_counts_by_zd[k][d] += count_sum;

		for (int m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
			i = words[m];
			w = d_state->m_words[i].m_word_index;
			m_word_counts_by_zw[k_old][w]--;
			m_word_counts_by_zw[k][w]++;
		}
		if (k == m_num_topics) // a new topic is created
				{
			m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_z.size() < m_num_topics + 1) {
				m_num_tables_by_z.push_back(0);
				m_word_counts_by_z.push_back(0);

				int* p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_word_counts_by_zd.push_back(p);

				p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
		}
	}
	delete[] counts;
	delete[] counts_copy;
}

void hdp_state::sample_word_assignment(doc_state* d_state, int i, bool remove,
		double_vec & q, double_vec & f) {
	if (remove)
		doc_state_update(d_state, i, -1);

	if ((int) q.size() < d_state->m_num_tables + 1)
		q.resize(2 * d_state->m_num_tables + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	int k, t, w;
	w = d_state->m_words[i].m_word_index;
	double f_new = m_gamma / m_size_vocab;
	for (k = 0; k < m_num_topics; k++) {
		f[k] = (m_word_counts_by_zw[k][w] + m_eta)
				/ (m_word_counts_by_z[k] + m_size_vocab * m_eta);
		f_new += m_num_tables_by_z[k] * f[k];
	}
	f_new = f_new / (m_total_num_tables + m_gamma);

	double total_q = 0.0, f_k = 0.0;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			k = d_state->m_table_to_topic[t];
			f_k = f[k];
		} else
			f_k = 0.0;

		total_q += d_state->m_word_counts_by_t[t] * f_k;
		q[t] = total_q;
	}
	total_q += m_alpha * f_new;
	q[d_state->m_num_tables] = total_q;

	double u = runiform() * total_q;
	for (t = 0; t < d_state->m_num_tables + 1; t++)
		if (u < q[t])
			break;

	d_state->m_words[i].m_table_assignment = t; // assign the new table

	if (t == d_state->m_num_tables) // this is a new table, we need get its k
			{
		if ((int) q.size() < m_num_topics + 1)
			q.resize(2 * m_num_topics + 1, 0.0);

		total_q = 0.0;
		for (k = 0; k < m_num_topics; k++) {
			total_q += m_num_tables_by_z[k] * f[k];
			q[k] = total_q;
		}
		total_q += m_gamma / m_size_vocab;
		q[m_num_topics] = total_q;
		u = runiform() * total_q;
		for (k = 0; k < m_num_topics + 1; k++)
			if (u < q[k])
				break;
		doc_state_update(d_state, i, +1, k);
	} else {
		doc_state_update(d_state, i, +1);
	}
}

void hdp_state::sample_word_assignment_online(doc_state* d_state, int i,
		bool remove, double_vec & q, double_vec & f) {
	if (remove)
		doc_state_update_online(d_state, i, -1);

	if ((int) q.size() < d_state->m_num_tables + 1)
		q.resize(2 * d_state->m_num_tables + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	int k, t, w;
	w = d_state->m_words[i].m_word_index;
	double f_new = m_gamma / m_size_vocab;
	for (k = 0; k < m_num_topics; k++) {
		f[k] = (m_word_counts_by_zw[k][w] + m_eta)
				/ (m_word_counts_by_z[k] + m_size_vocab * m_eta);
		f_new += m_num_tables_by_z[k] * f[k];
	}
	f_new = f_new / (m_total_num_tables + m_gamma);

	double total_q = 0.0, f_k = 0.0;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			k = d_state->m_table_to_topic[t];
			f_k = f[k];
		} else
			f_k = 0.0;

		total_q += d_state->m_word_counts_by_t[t] * f_k;
		q[t] = total_q;
	}
	total_q += m_alpha * f_new;
	q[d_state->m_num_tables] = total_q;

	double u = runiform() * total_q;
	for (t = 0; t < d_state->m_num_tables + 1; t++)
		if (u < q[t])
			break;

	d_state->m_words[i].m_table_assignment = t; // assign the new table

	if (t == d_state->m_num_tables) // this is a new table, we need get its k
			{
		if ((int) q.size() < m_num_topics + 1)
			q.resize(2 * m_num_topics + 1, 0.0);

		total_q = 0.0;
		for (k = 0; k < m_num_topics; k++) {
			total_q += m_num_tables_by_z[k] * f[k];
			q[k] = total_q;
		}
		total_q += m_gamma / m_size_vocab;
		q[m_num_topics] = total_q;
		u = runiform() * total_q;
		for (k = 0; k < m_num_topics + 1; k++)
			if (u < q[k])
				break;
		doc_state_update_online(d_state, i, +1, k);
	} else {
		doc_state_update_online(d_state, i, +1);
	}
}

// k is only provided when m_table_to_topic doesn't have that
void hdp_state::doc_state_update(doc_state* d_state, int i, int update, int k) {
	int d, w, t;
	d = d_state->m_doc_id;
	w = d_state->m_words[i].m_word_index;
	//k = d_state->m_words[i].m_topic_assignment;
	t = d_state->m_words[i].m_table_assignment;
	if (k < 0)
		k = d_state->m_table_to_topic[t];
	assert(k >= 0);

	d_state->m_word_counts_by_t[t] += update;

	m_word_counts_by_z[k] += update;
	m_word_counts_by_zw[k][w] += update;
	m_word_counts_by_zd[k][d] += update;

	if (update == -1 && d_state->m_word_counts_by_t[t] == 0) /// this table becomes empty
			{
		m_total_num_tables--;
		m_num_tables_by_z[k]--;
		d_state->m_table_to_topic[t] = -1;
		/// m_num_topics, no need to change at this moment
	}

	if (update == 1 && d_state->m_word_counts_by_t[t] == 1) /// a new table is created
			{
		if (t == d_state->m_num_tables)
			d_state->m_num_tables++; // create a new table
		d_state->m_table_to_topic[t] = k; // mapping the table

		m_num_tables_by_z[k]++;          // adding the table to mixture k
		m_total_num_tables++;

		if ((int) d_state->m_table_to_topic.size()
				< d_state->m_num_tables + 1) {
			d_state->m_table_to_topic.push_back(-1);
			d_state->m_word_counts_by_t.push_back(0);
		}
		if (k == m_num_topics) // used to k == m_num_topics
				{
			assert(m_word_counts_by_z[k] == 1);
			if (k == m_num_topics)
				m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_z.size() < m_num_topics + 1) {
				m_num_tables_by_z.push_back(0);
				m_word_counts_by_z.push_back(0);

				int* p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_word_counts_by_zd.push_back(p);

				p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
		}
	}
}

// k is only provided when m_table_to_topic doesn't have that
void hdp_state::doc_state_update_online(doc_state* d_state, int i, int update,
		int k) {
	int d, w, t;
	d = d_state->m_doc_id;
	w = d_state->m_words[i].m_word_index;
	//k = d_state->m_words[i].m_topic_assignment;
	t = d_state->m_words[i].m_table_assignment;
	if (k < 0)
		k = d_state->m_table_to_topic[t];
	assert(k >= 0);

	if (update == -1 && d_state->m_word_counts_by_t[t] == 0) /// this table becomes empty
			{
		m_total_num_tables--;
		m_num_tables_by_z[k]--;
		d_state->m_table_to_topic[t] = -1;

		/// m_num_topics, no need to change at this moment
	}

	d_state->m_word_counts_by_t[t] += update;

	if (update == 1 && d_state->m_word_counts_by_t[t] == 1) /// a new table is created
			{
		if (k == m_num_topics) // used to k == m_num_topics
				{
			//assert(m_word_counts_by_z[k] == 1);
			if (k == m_num_topics)
				m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_z.size() < m_num_topics + 1) {
				m_num_tables_by_z.push_back(0);
				m_word_counts_by_z.push_back(0);

				//int* p = new int [m_num_docs];
				//memset(p, 0, sizeof(int)*m_num_docs);
				//m_word_counts_by_zd.push_back(p);

				int * p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
			if ((int) m_word_counts_by_zd_online.size() < m_num_topics) {
				m_word_counts_by_zd_online.push_back(0);
			}
		}

		if (t == d_state->m_num_tables)
			d_state->m_num_tables++; // create a new table
		d_state->m_table_to_topic[t] = k; // mapping the table

		m_num_tables_by_z[k]++;          // adding the table to mixture k
		m_total_num_tables++;

		if ((int) d_state->m_table_to_topic.size()
				< d_state->m_num_tables + 1) {
			d_state->m_table_to_topic.push_back(-1);
			d_state->m_word_counts_by_t.push_back(0);
		}
	}

	m_word_counts_by_z[k] += update;
	m_word_counts_by_zw[k][w] += update;
	m_word_counts_by_zd_online[k] += update;
}

void hdp_state::save_state(char * name) {
	char filename[500];

	// save the topic words counts
	sprintf(filename, "%s-topics.dat", name);
	FILE* file = fopen(filename, "w");

	for (int k = 0; k < m_num_topics; k++) {
		for (int w = 0; w < m_size_vocab; w++)
			fprintf(file, "%05d ", m_word_counts_by_zw[k][w]);
		fprintf(file, "\n");
	}
	fclose(file);

	sprintf(filename, "%s-word-assignments.dat", name);
	file = fopen(filename, "w");
	fprintf(file, "d w z t\n");
	int w, k, t;
	for (int d = 0; d < m_num_docs; d++) {
		doc_state* d_state = m_doc_states[d];
		int doc_id = d_state->m_doc_id;
		for (int i = 0; i < d_state->m_doc_length; i++) {
			w = d_state->m_words[i].m_word_index;
			t = d_state->m_words[i].m_table_assignment;
			k = d_state->m_table_to_topic[t];
			fprintf(file, "%d %d %d %d\n", doc_id, w, k, t);
		}
	}
	fclose(file);
}

void hdp_state::save_state_ex(char * name, int sample_num) {

	if (sample_num == -1) {
		save_state_to_file(name);
	} else {
		char * file_name = new char[strlen(name) + 5 * sizeof(char)];
		strcpy(file_name, name);
		char* postfix = new char[5];
		postfix[0] = '_';
		postfix[1] = '0' + sample_num / 100;
		postfix[2] = '0' + (sample_num / 10) % 10;
		postfix[3] = '0' + sample_num % 10;
		postfix[4] = '\0';

		file_name = strcat(file_name, postfix);
		save_state_to_file(file_name);
		delete[] file_name;
		delete[] postfix;
	}

}

void hdp_state::save_state_to_file(char * name) {
	FILE * file = fopen(name, "wb");
	fwrite(&m_size_vocab, sizeof(int), 1, file);
	fwrite(&m_total_words, sizeof(int), 1, file);
	fwrite(&m_num_topics, sizeof(int), 1, file);
	fwrite(&m_total_num_tables, sizeof(int), 1, file);

	fwrite(&m_eta, sizeof(double), 1, file);
	fwrite(&m_gamma, sizeof(double), 1, file);
	fwrite(&m_alpha, sizeof(double), 1, file);

	for (int k = 0; k < m_num_topics; k++) {
		fwrite(&(m_num_tables_by_z[k]), sizeof(int), 1, file);
		fwrite(&(m_word_counts_by_z[k]), sizeof(int), 1, file);
		fwrite(m_word_counts_by_zw[k], sizeof(int), m_size_vocab, file);
	}
	fclose(file);
}

void hdp_state::load_state_ex(char * name, int sample_num) {
	if (sample_num == -1) {
		load_state_from_file(name);
	} else {
		char * file_name = new char[strlen(name) + 5 * sizeof(char)];
		strcpy(file_name, name);
		char* postfix = new char[5];
		postfix[0] = '_';
		postfix[1] = '0' + sample_num / 100;
		postfix[2] = '0' + (sample_num / 10) % 10;
		postfix[3] = '0' + sample_num % 10;
		postfix[4] = '\0';

		file_name = strcat(file_name, postfix);
		load_state_from_file(file_name);
		delete[] file_name;
		delete[] postfix;
	}
}

void hdp_state::load_state_from_file(char * name) {
	free_state();
	FILE * file = fopen(name, "rb");
	fread(&m_size_vocab, sizeof(int), 1, file);
	fread(&m_total_words, sizeof(int), 1, file);
	fread(&m_num_topics, sizeof(int), 1, file);
	fread(&m_total_num_tables, sizeof(int), 1, file);

	fread(&m_eta, sizeof(double), 1, file);
	fread(&m_gamma, sizeof(double), 1, file);
	fread(&m_alpha, sizeof(double), 1, file);

	m_num_tables_by_z.resize(m_num_topics);
	m_word_counts_by_z.resize(m_num_topics);
	m_word_counts_by_zw.resize(m_num_topics);
	for (int k = 0; k < m_num_topics; k++) {
		fread(&(m_num_tables_by_z[k]), sizeof(int), 1, file);
		fread(&(m_word_counts_by_z[k]), sizeof(int), 1, file);

		m_word_counts_by_zw[k] = new int[m_size_vocab];
		fread(m_word_counts_by_zw[k], sizeof(int), m_size_vocab, file);
	}
	fclose(file);

}

hdp_state_dynamic::hdp_state_dynamic() {
	m_doc_states = NULL;
	m_size_vocab = 0;
	m_total_words = 0;
	;
	m_num_docs = 0;
	m_num_topics = 0;
	m_total_num_tables = 0;

	m_word_counts_by_z.clear();
	m_num_tables_by_z.clear();
	m_word_counts_by_zd.clear();
	m_num_tables_by_zd.clear();
	m_total_num_tables_by_d.clear();
	m_word_counts_by_zw.clear();

	m_word_counts_by_zd_online.clear();
	m_num_tables_by_zd_online.clear();
}

hdp_state_dynamic::~hdp_state_dynamic() {
	free_state();
}

void hdp_state_dynamic::free_state() {
	if (m_doc_states != NULL) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
		}
		delete[] m_doc_states;
	}
	m_doc_states = NULL;

	m_size_vocab = 0;
	m_total_words = 0;
	m_num_docs = 0;
	m_num_topics = 0;
	m_num_topics_before_update = 0;
	m_total_num_tables_by_d_online = 0;
	m_total_num_tables_by_d_for_prev_doc_online = 0;
	m_total_num_tables = 0;

	m_word_counts_by_z.clear();
	m_num_tables_by_z.clear();
	m_total_num_tables_by_d.clear();

	free_vec_ptr(m_word_counts_by_zd);
	free_vec_ptr(m_num_tables_by_zd);
	free_vec_ptr(m_word_counts_by_zw);

	m_word_counts_by_zd_online.clear();
	m_num_tables_by_zd_online.clear();
	m_num_tables_by_zd_for_prev_doc_online.clear();
}

void hdp_state_dynamic::set_vocab_size(int _vocab_size) {
	m_size_vocab = _vocab_size;
}

void hdp_state_dynamic::setup_state_from_corpus(const corpus * c) {
	m_size_vocab = max(c->size_vocab, m_size_vocab);
	m_total_words += c->total_words;
	m_num_docs = c->num_docs;
	m_doc_states = new doc_state *[m_num_docs];

	for (unsigned int d = 0; d < c->docs.size(); d++) {
		document * doc = c->docs[d];
		doc_state * d_state = new doc_state();
		m_doc_states[d] = d_state;
		d_state->setup_state_from_doc(doc);
	}
}

void hdp_state_dynamic::setup_doc_info_from_document(const document * doc) {
	int proposed_vocab_size = doc->get_max_word_id() + 1;
	m_size_vocab = max(m_size_vocab, proposed_vocab_size);
	m_total_words += doc->total;

	if (!(m_doc_states == NULL)) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
			d_state = NULL;
		}
		delete[] m_doc_states;
		m_doc_states = NULL;
	}
	m_num_docs = 1;
	m_doc_states = new doc_state *[1];

	doc_state * d_state = new doc_state();
	m_doc_states[0] = d_state;
	d_state->setup_state_from_doc(doc);
}

void hdp_state_dynamic::allocate_initial_space() {
	// training
	if (m_word_counts_by_z.size() == 0) {
		m_total_num_tables_by_d.resize(m_num_docs, 0);

		m_word_counts_by_z.resize(INIT_SIZE, 0);
		m_num_tables_by_z.resize(INIT_SIZE, 0);
		m_word_counts_by_zd.resize(INIT_SIZE, NULL);
		m_num_tables_by_zd.resize(INIT_SIZE, NULL);
		m_word_counts_by_zw.resize(INIT_SIZE, NULL);
		int * p = NULL;
		for (int k = 0; k < INIT_SIZE; k++) {
			p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_word_counts_by_zd[k] = p;

			p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_num_tables_by_zd[k] = p;

			p = new int[m_size_vocab];
			memset(p, 0, sizeof(int) * m_size_vocab);
			m_word_counts_by_zw[k] = p;
		}
	} else // testing
	{
		while ((int) m_word_counts_by_z.size() < m_num_topics + 1) {
			m_word_counts_by_z.push_back(0);
			m_num_tables_by_z.push_back(0);

			int* p = new int[m_size_vocab];
			memset(p, 0, sizeof(int) * m_size_vocab);
			m_word_counts_by_zw.push_back(p);
		}
		while ((int) m_word_counts_by_zd.size() < m_num_topics + 1) {
			int* p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_word_counts_by_zd.push_back(p);
		}
		while ((int) m_num_tables_by_zd.size() < m_num_topics + 1) {
			int* p = new int[m_num_docs];
			memset(p, 0, sizeof(int) * m_num_docs);
			m_num_tables_by_zd.push_back(p);
		}
		m_total_num_tables_by_d.resize(m_num_docs, 0);
	}
}

void hdp_state_dynamic::allocate_space_for_online_updates() {
	m_word_counts_by_zd_online.resize(m_num_topics, 0);
	m_num_tables_by_zd_online.resize(m_num_topics, 0);

	m_total_num_tables_by_d_online = 0;

	m_num_topics_before_update = m_num_topics;
}

void hdp_state_dynamic::free_doc_info() {
	if (m_doc_states != NULL) {
		for (int d = 0; d < m_num_docs; d++) {
			doc_state * d_state = m_doc_states[d];
			delete d_state;
		}
		delete[] m_doc_states;
	}
	m_doc_states = NULL;

	m_num_docs = 0;

	for (int k = 0; k < m_word_counts_by_zd.size(); ++k) {
		delete[] m_word_counts_by_zd[k];
	}
	m_word_counts_by_zd.clear();

	for (int k = 0; k < m_num_tables_by_zd.size(); ++k) {
		delete[] m_num_tables_by_zd[k];
	}
	m_num_tables_by_zd.clear();

	m_total_num_tables_by_d.clear();
}

void hdp_state_dynamic::free_counts_update() {
	m_word_counts_by_zd_online.clear();
	m_num_tables_by_zd_online.clear();
	m_num_tables_by_zd_for_prev_doc_online.clear();

	m_total_num_tables_by_d_online = 0;
	m_total_num_tables_by_d_for_prev_doc_online = 0;
	m_num_topics_before_update = 0;
}

void hdp_state_dynamic::free_online_info() {
	free_doc_info();
	free_counts_update();
}

void hdp_state_dynamic::iterate_gibbs_state(bool not_init,
		hdp_hyperparameter * hdp_hyperparam, bool table_sampling) {
	double_vec q;
	double_vec f;
	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		reset_counts_up_to_date(j);
		d_state = m_doc_states[j];
		for (int i = 0; i < d_state->m_doc_length; i++) {
			sample_word_assignment(d_state, i, not_init, q, f);
		}
		if (table_sampling)
			sample_tables(d_state, not_init, q, f);

	}
	compact_hdp_state();
}

void hdp_state_dynamic::iterate_gibbs_state_online(bool remove,
		hdp_hyperparameter* hdp_hyperparam, bool table_sampling) {
	double_vec q;
	double_vec f;
	doc_state* d_state = m_doc_states[0];

	for (int i = 0; i < d_state->m_doc_length; i++) {
		sample_word_assignment_online(d_state, i, remove, q, f);
	}
	if (table_sampling)
		sample_tables_online(d_state, q, f);

	compact_hdp_state_online();
}

void hdp_state_dynamic::reset_counts_up_to_date(int doc_id) {
	if (m_total_num_tables_by_d.size() > 0) {
		if (doc_id == 0) {
			// reset count for the first document
			m_total_num_tables = m_total_num_tables_by_d[doc_id];
			m_num_tables_by_z.resize(m_num_tables_by_zd.size(), 0);
			for (int k = 0; k < m_num_tables_by_zd.size(); ++k) {
				assert(m_num_tables_by_zd[k] != NULL);
				m_num_tables_by_z[k] = m_num_tables_by_zd[k][doc_id];
			}
		} else {
			// update count
			m_total_num_tables += m_total_num_tables_by_d[doc_id];
			for (int k = 0; k < m_num_tables_by_z.size(); ++k) {
				assert(m_num_tables_by_zd[k] != NULL);
				m_num_tables_by_z[k] += m_num_tables_by_zd[k][doc_id];
			}
			for (int k = m_num_tables_by_z.size();
					k < m_num_tables_by_zd.size(); ++k) {
				assert(m_num_tables_by_zd[k] != NULL);
				m_num_tables_by_z.push_back(m_num_tables_by_zd[k][doc_id]);
			}
		}
	} else {
		m_total_num_tables = 0;
		m_num_tables_by_z.resize(INIT_SIZE, 0);
	}
}

void hdp_state_dynamic::sample_tables(doc_state* d_state, bool not_init,
		double_vec & q, double_vec & f) {
	vector<int*> words_by_t;
	words_by_t.resize(d_state->m_num_tables, NULL);
	int* p = NULL;
	int t, word;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			p = new int[d_state->m_word_counts_by_t[t]];
			words_by_t[t] = p;
		}
	}

	int* i_by_t = new int[d_state->m_num_tables];
	memset(i_by_t, 0, sizeof(int) * d_state->m_num_tables);

	for (int i = 0; i < d_state->m_doc_length; i++) {
		word = d_state->m_words[i].m_word_index;
		t = d_state->m_words[i].m_table_assignment;
		words_by_t[t][i_by_t[t]] = i; // save the id for later use
		i_by_t[t]++;
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0)
			sample_table_assignment(d_state, t, not_init, words_by_t[t], q, f);
		//eles no needs, since there is no data there
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		p = words_by_t[t];
		delete[] p;
	}
	delete[] i_by_t;
}

void hdp_state_dynamic::sample_tables_online(doc_state* d_state, double_vec & q,
		double_vec & f) {
	vector<int*> words_by_t;
	words_by_t.resize(d_state->m_num_tables, NULL);
	int* p = NULL;
	int t, word;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			p = new int[d_state->m_word_counts_by_t[t]];
			words_by_t[t] = p;
		}
	}

	int* i_by_t = new int[d_state->m_num_tables];
	memset(i_by_t, 0, sizeof(int) * d_state->m_num_tables);

	for (int i = 0; i < d_state->m_doc_length; i++) {
		word = d_state->m_words[i].m_word_index;
		t = d_state->m_words[i].m_table_assignment;
		words_by_t[t][i_by_t[t]] = i; // save the id for later use
		i_by_t[t]++;
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0)
			sample_table_assignment_online(d_state, t, words_by_t[t], q, f);
		//eles no needs, since there is no data there
	}

	for (t = 0; t < d_state->m_num_tables; t++) {
		p = words_by_t[t];
		delete[] p;
	}
	delete[] i_by_t;
}

void hdp_state_dynamic::sample_table_assignment(doc_state* d_state, int t,
		bool not_init, int * words, double_vec & q, double_vec & f) {
	//number of tables won't change at all
	int i, w, k, m, k_old, d;
	int* counts = new int[m_size_vocab];
	int* counts_copy = new int[m_size_vocab];
	memset(counts_copy, 0, sizeof(int) * m_size_vocab);
	int count_sum = d_state->m_word_counts_by_t[t];

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		counts_copy[w]++;
	}
	memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);

	// compute the the log prob of being at a new cluster
	double f_new = lgamma(m_size_vocab * m_eta)
			- lgamma(count_sum + m_size_vocab * m_eta);

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		if (counts[w] > 0) {
			f_new += lgamma(counts[w] + m_eta) - lgamma(m_eta);
			counts[w] = 0;
		}
	}

	if ((int) q.size() < m_num_topics + 1)
		q.resize(2 * m_num_topics + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	double p_k_new = 0;

	k_old = d_state->m_table_to_topic[t];

	vector<double> topic_prior_probability;
	topic_prior_probability.assign(m_num_topics + 1, 0);

	vector<double> next_topics_probability;
	next_topics_probability.assign(m_num_topics + 1, 0);

	double new_topic_prior = compute_prior_probability_for_topic_sampling(
			d_state->m_doc_id, m_num_topics, false);
	topic_prior_probability[m_num_topics] = new_topic_prior;

	double next_topic_new = 0;
	if ((not_init) && (topic_prior_probability[m_num_topics] > EPS)) {
		next_topic_new = compute_next_table_assignment_likelihood(
				d_state->m_doc_id, k_old, -1, false);
	}
	next_topics_probability[m_num_topics] = next_topic_new;

	for (k = 0; k < m_num_topics; k++) {
		if (k == k_old) {
			f[k] = lgamma(
					m_size_vocab * m_eta + m_word_counts_by_z[k] - count_sum)
					- lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k]);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(m_eta + m_word_counts_by_zw[k][w])
							- lgamma(
									m_eta + m_word_counts_by_zw[k][w]
											- counts[w]);
					counts[w] = 0;
				}
			}

			double topic_prior = compute_prior_probability_for_topic_sampling(
					d_state->m_doc_id, k, true);
			topic_prior_probability[k] = topic_prior;

			double next_topic = 0;
			if ((not_init) && (topic_prior_probability[k] > EPS)) {
				next_topic = compute_next_table_assignment_likelihood(
						d_state->m_doc_id, k_old, k, false);
			}
			next_topics_probability[k] = next_topic;
		} else {
			f[k] = lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k])
					- lgamma(
							m_size_vocab * m_eta + m_word_counts_by_z[k]
									+ count_sum);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(
							m_eta + m_word_counts_by_zw[k][w] + counts[w])
							- lgamma(m_eta + m_word_counts_by_zw[k][w]);
					counts[w] = 0;
				}
			}

			double topic_prior = compute_prior_probability_for_topic_sampling(
					d_state->m_doc_id, k, false);
			topic_prior_probability[k] = topic_prior;

			double next_topic = 0;
			if ((not_init) && (topic_prior_probability[k] > EPS)) {
				next_topic = compute_next_table_assignment_likelihood(
						d_state->m_doc_id, k_old, k, false);
			}
			next_topics_probability[k] = next_topic;
		}
	}
	if (not_init) {
		log_normalize_non_zero_values(next_topics_probability,
				m_num_topics + 1);
	}
	for (k = 0; k < m_num_topics; ++k) {
		if ((k != k_old) || (topic_prior_probability[k] > EPS)) {
			q[k] = log(topic_prior_probability[k]) + next_topics_probability[k]
					+ f[k];
		} else {
			q[k] = INF;
		}
	}
	q[m_num_topics] = log(topic_prior_probability[m_num_topics])
			+ next_topics_probability[m_num_topics] + f_new;

	//normalizing in log space for sampling
	log_normalize(q, m_num_topics + 1);
	q[0] = exp(q[0]);
	double total_q = q[0];
	for (k = 1; k < m_num_topics + 1; k++) {
		total_q += exp(q[k]);
		q[k] = total_q;
	}

	double u = runiform() * total_q;
	for (k = 0; k < m_num_topics + 1; k++) {
		if (u < q[k])
			break;
	}

	if (k != k_old) // status doesn't change, but k could change
			{
		d = d_state->m_doc_id;

		/// reassign the topic to current table
		d_state->m_table_to_topic[t] = k;

		/// update the statistics by removing the table t from topic k_old
		m_num_tables_by_zd[k_old][d]--;
		m_num_tables_by_z[k_old]--;
		m_word_counts_by_z[k_old] -= count_sum;
		m_word_counts_by_zd[k_old][d] -= count_sum;

		if (k == m_num_topics) // a new topic is created
				{
			m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_zd.size() < m_num_topics + 1) {
				m_word_counts_by_z.push_back(0);
				m_num_tables_by_z.push_back(0);

				int* p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_word_counts_by_zd.push_back(p);

				p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_num_tables_by_zd.push_back(p);

				p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
		}

		/// update the statistics by adding the table t to topic k
		m_num_tables_by_zd[k][d]++;
		m_num_tables_by_z[k]++;
		m_word_counts_by_z[k] += count_sum;
		m_word_counts_by_zd[k][d] += count_sum;

		for (int m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
			i = words[m];
			w = d_state->m_words[i].m_word_index;
			m_word_counts_by_zw[k_old][w]--;
			m_word_counts_by_zw[k][w]++;
		}
	}
	delete[] counts;
	delete[] counts_copy;
}

double hdp_state_dynamic::compute_next_doc_table_assignment_likelihood(
		int doc_id, int old_topic_id, int new_topic_id, bool new_table) const {
	if (doc_id == m_num_docs - 1) {
		// there is no next document, its likelihood = 1 => log likelihood = 0;
		return 0;
	}

	double likelihood =
			compute_next_table_assignment_likelihood_for_non_born_topics(doc_id,
					m_num_tables_by_z, old_topic_id, new_topic_id, true);

	return likelihood;
}

double hdp_state_dynamic::compute_next_table_assignment_likelihood_for_non_born_topics(
		int doc_id, const vector<int>& updating_num_table_by_z,
		int old_topic_id, int new_topic_id,
		bool update_counts_for_cur_doc) const {
	double likelihood = 0;

	if ((old_topic_id == new_topic_id)
			|| ((old_topic_id == -1) && (new_topic_id == -1))) {
		for (int k = 0; k < m_num_topics; ++k) {
			if (updating_num_table_by_z[k] > 0) {
				for (int n = 0; n < m_num_tables_by_zd[k][doc_id + 1]; ++n) {
					likelihood += log(
							m_num_tables_by_zd[k][doc_id] + n
									+ m_delta
											* (updating_num_table_by_z[k] + n));
				}
			}
		}
	} else {
		for (int k = 0; k < m_num_topics; ++k) {
			if (k == old_topic_id) {
				if (updating_num_table_by_z[k] - 1 > 0) {
					int m_tables_in_current_doc = m_num_tables_by_zd[k][doc_id];
					if (update_counts_for_cur_doc) {
						--m_tables_in_current_doc;
					}
					for (int n = 0; n < m_num_tables_by_zd[k][doc_id + 1];
							++n) {
						likelihood += log(
								m_tables_in_current_doc + n
										+ m_delta
												* (updating_num_table_by_z[k]
														- 1 + n));
					}
				}
			} else if (k == new_topic_id) {
				if (updating_num_table_by_z[k] + 1 > 0) {
					int m_tables_in_current_doc = m_num_tables_by_zd[k][doc_id];
					if (update_counts_for_cur_doc) {
						++m_tables_in_current_doc;
					}
					for (int n = 0; n < m_num_tables_by_zd[k][doc_id + 1];
							++n) {
						likelihood += log(
								m_tables_in_current_doc + n
										+ m_delta
												* (updating_num_table_by_z[k]
														+ 1 + n));
					}
				}
			} else {
				if (updating_num_table_by_z[k] > 0) {
					for (int n = 0; n < m_num_tables_by_zd[k][doc_id + 1];
							++n) {
						likelihood += log(
								m_num_tables_by_zd[k][doc_id] + n
										+ m_delta
												* (updating_num_table_by_z[k]
														+ n));
					}
				}
			}
		}
	}

	return likelihood;
}

double hdp_state_dynamic::compute_next_table_assignment_likelihood(int doc_id,
		int old_topic_id, int new_topic_id, bool new_table) const {
	if (doc_id == m_num_docs - 1) {
		return 0; // log(0) = 1
	}
	double next_doc_probability = compute_next_doc_table_assignment_likelihood(
			doc_id, old_topic_id, new_topic_id, new_table);

	double other_next_docs_probability = 0;
	double result = next_doc_probability + other_next_docs_probability;

	return result;
}

void hdp_state_dynamic::sample_table_assignment_online(doc_state* d_state,
		int t, int* words, double_vec & q, double_vec & f) {
	//number of tables won't change at all
	int i, w, k, m, k_old, d;
	int* counts = new int[m_size_vocab];
	int* counts_copy = new int[m_size_vocab];
	memset(counts_copy, 0, sizeof(int) * m_size_vocab);
	int count_sum = d_state->m_word_counts_by_t[t];

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		counts_copy[w]++;
	}
	memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);

	// compute the the log prob of being at a new cluster
	double f_new = lgamma(m_size_vocab * m_eta)
			- lgamma(count_sum + m_size_vocab * m_eta);

	for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
		i = words[m];
		w = d_state->m_words[i].m_word_index;
		if (counts[w] > 0) {
			f_new += lgamma(counts[w] + m_eta) - lgamma(m_eta);
			counts[w] = 0;
		}
	}

	if ((int) q.size() < m_num_topics + 1)
		q.resize(2 * m_num_topics + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	double topic_prior_new =
			compute_prior_probability_for_topic_sampling_online(m_num_topics,
					false);
	q[m_num_topics] = log(topic_prior_new) + f_new;

	k_old = d_state->m_table_to_topic[t];

	for (k = 0; k < m_num_topics; k++) {
		if (k == k_old) {
			f[k] = lgamma(
					m_size_vocab * m_eta + m_word_counts_by_z[k] - count_sum)
					- lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k]);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(m_eta + m_word_counts_by_zw[k][w])
							- lgamma(
									m_eta + m_word_counts_by_zw[k][w]
											- counts[w]);
					counts[w] = 0;
				}
			}

			double topic_prior =
					compute_prior_probability_for_topic_sampling_online(k,
							true);

			if (topic_prior < EPS) {
				q[k] = INF; // make it extremely small as log(0)
			} else {
				q[k] = log(topic_prior) + f[k];
			}
		} else {
			f[k] = lgamma(m_size_vocab * m_eta + m_word_counts_by_z[k])
					- lgamma(
							m_size_vocab * m_eta + m_word_counts_by_z[k]
									+ count_sum);

			memcpy(counts, counts_copy, sizeof(int) * m_size_vocab);
			for (m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
				i = words[m];
				w = d_state->m_words[i].m_word_index;
				if (counts[w] > 0) {
					f[k] += lgamma(
							m_eta + m_word_counts_by_zw[k][w] + counts[w])
							- lgamma(m_eta + m_word_counts_by_zw[k][w]);
					counts[w] = 0;
				}
			}

			double topic_prior =
					compute_prior_probability_for_topic_sampling_online(k,
							false);

			q[k] = log(topic_prior) + f[k];
		}
	}
	//normalizing in log space for sampling
	log_normalize(q, m_num_topics + 1);
	q[0] = exp(q[0]);
	double total_q = q[0];
	for (k = 1; k < m_num_topics + 1; k++) {
		total_q += exp(q[k]);
		q[k] = total_q;
	}

	double u = runiform() * total_q;
	for (k = 0; k < m_num_topics + 1; k++)
		if (u < q[k])
			break;

	if (k != k_old) // status doesn't change, but k could change
			{
		d = d_state->m_doc_id;

		/// reassign the topic to current table
		d_state->m_table_to_topic[t] = k;

		/// update the statistics by removing the table t from topic k_old
		m_num_tables_by_zd_online[k_old]--;
		m_word_counts_by_z[k_old] -= count_sum;
		m_num_tables_by_z[k_old]--;
		m_word_counts_by_zd_online[k_old] -= count_sum;

		if (k == m_num_topics) // a new topic is created
				{
			m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_zd_online.size() < m_num_topics + 1) {
				m_num_tables_by_zd_online.push_back(0);
				m_word_counts_by_z.push_back(0);
				m_num_tables_by_z.push_back(0);

				int * p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
			if ((int) m_word_counts_by_zd_online.size() < m_num_topics) {
				m_word_counts_by_zd_online.push_back(0);
			}
		}

		/// update the statistics by adding the table t to topic k
		m_num_tables_by_zd_online[k]++;
		m_num_tables_by_z[k]++;
		m_word_counts_by_z[k] += count_sum;
		m_word_counts_by_zd_online[k] += count_sum;

		for (int m = 0; m < d_state->m_word_counts_by_t[t]; m++) {
			i = words[m];
			w = d_state->m_words[i].m_word_index;
			m_word_counts_by_zw[k_old][w]--;
			m_word_counts_by_zw[k][w]++;
		}
	}
	delete[] counts;
	delete[] counts_copy;
}

void hdp_state_dynamic::sample_word_assignment(doc_state* d_state, int i,
		bool not_init, double_vec & q, double_vec & f) {
	if (not_init)
		doc_state_update(d_state, i, -1);

	if ((int) q.size() < d_state->m_num_tables + 1)
		q.resize(2 * d_state->m_num_tables + 1, 0.0);

	if ((int) f.size() <= m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	vector<double> topic_sampling_probability;
	topic_sampling_probability.assign(m_num_topics + 1, 0);

	vector<double> prior_topic_probability;
	prior_topic_probability.assign(m_num_topics + 1, 0);

	vector<double> next_doc_topic_assignment_probability;
	next_doc_topic_assignment_probability.assign(m_num_topics + 1, 0);

	int k, t, w;
	w = d_state->m_words[i].m_word_index;

	double topic_prior_normalisation_constant = 0;

	double f_new = 0;
	for (k = 0; k < m_num_topics; k++) {
		f[k] = (m_word_counts_by_zw[k][w] + m_eta)
				/ (m_word_counts_by_z[k] + m_size_vocab * m_eta);

		double topic_prior = compute_prior_probability_for_topic_sampling(
				d_state->m_doc_id, k, false);
		prior_topic_probability[k] = topic_prior;

		double next_topics_probability = 0;
		if ((not_init) && (prior_topic_probability[k] > EPS)) {
			next_topics_probability = compute_next_table_assignment_likelihood(
					d_state->m_doc_id, -1, k, true);
		}
		next_doc_topic_assignment_probability[k] = next_topics_probability;
	}

	double topic_new_prior = compute_prior_probability_for_topic_sampling(
			d_state->m_doc_id, m_num_topics, false);
	prior_topic_probability[m_num_topics] = topic_new_prior;

	double next_topics_probability_new = 0;
	if ((not_init) && (prior_topic_probability[m_num_topics] > EPS)) {
		next_topics_probability_new = compute_next_table_assignment_likelihood(
				d_state->m_doc_id, -1, -1, true);
	}
	next_doc_topic_assignment_probability[m_num_topics] =
			next_topics_probability_new;

	if (not_init) {
		log_normalize_non_zero_values(next_doc_topic_assignment_probability,
				m_num_topics + 1);
	}

	for (k = 0; k < m_num_topics; ++k) {
		double topic_probability_term = prior_topic_probability[k]
				* exp(next_doc_topic_assignment_probability[k]);
		topic_prior_normalisation_constant += topic_probability_term;
		f_new += topic_probability_term * f[k];
		topic_sampling_probability[k] = f_new;
	}

	double new_topic_probability_term = prior_topic_probability[m_num_topics]
			* exp(next_doc_topic_assignment_probability[m_num_topics]);
	topic_prior_normalisation_constant += new_topic_probability_term;
	f_new += new_topic_probability_term * 1.0 / m_size_vocab;
	topic_sampling_probability[m_num_topics] = f_new;

	double f_new_normalised = f_new / topic_prior_normalisation_constant;

	double total_q = 0.0, f_k = 0.0;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			k = d_state->m_table_to_topic[t];
			f_k = f[k];
		} else
			f_k = 0.0;

		total_q += d_state->m_word_counts_by_t[t] * f_k;
		q[t] = total_q;
	}
	total_q += m_alpha * f_new_normalised;
	q[d_state->m_num_tables] = total_q;

	double u = runiform() * total_q;
	for (t = 0; t < d_state->m_num_tables + 1; t++) {
		if (u < q[t]) {
			break;
		}
	}

	d_state->m_words[i].m_table_assignment = t; // assign the new table

	if (t == d_state->m_num_tables) // this is a new table, we need get its k
			{
		double u = runiform() * f_new;
		int new_topic = -1;
		for (int l = 0; l < m_num_topics + 1; l++) {
			if (u < topic_sampling_probability[l]) {
				new_topic = l;
				break;
			}
		}

		doc_state_update(d_state, i, +1, new_topic);
	} else {
		doc_state_update(d_state, i, +1);
	}
}

int hdp_state_dynamic::sample_new_topic_for_new_table_online(
		const vector<double>& likelihood_per_topic_assignment) const {
	vector<double> q;
	q.assign(m_num_topics + 1, 0);

	double total_q = 0.0;
	for (int k = 0; k < m_num_topics; k++) {
		double prior_probability =
				compute_prior_probability_for_topic_sampling_online(k, false);

		total_q += prior_probability * likelihood_per_topic_assignment[k];
		q[k] = total_q;
	}

	double prior_probability =
			compute_prior_probability_for_topic_sampling_online(m_num_topics,
					false);

	total_q += prior_probability
			* likelihood_per_topic_assignment[m_num_topics];
	q[m_num_topics] = total_q;
	double u = runiform() * total_q;
	int new_topic = -1;
	for (int k = 0; k < m_num_topics + 1; k++) {
		if (u < q[k]) {
			new_topic = k;
			break;
		}
	}

	return new_topic;
}

double hdp_state_dynamic::compute_prior_probability_for_topic_sampling(
		int doc_id, int topic_id, bool decrease_table_counts) const {
	assert(topic_id <= m_num_topics);
	if (topic_id == m_num_topics) {
		return m_gamma;
	}
	int num_tables_from_prev_doc = 0;
	if (doc_id > 0) {
		num_tables_from_prev_doc = m_num_tables_by_zd[topic_id][doc_id - 1];
	}
	int m_num_tables_for_cur_doc = m_num_tables_by_zd[topic_id][doc_id];
	int m_total_num_tables_for_topic = m_num_tables_by_z[topic_id];
	if (decrease_table_counts) {
		--m_num_tables_for_cur_doc;
		--m_total_num_tables_for_topic;
	}
	double prior_for_topic = m_num_tables_for_cur_doc + num_tables_from_prev_doc
			+ m_delta * m_total_num_tables_for_topic;

	return prior_for_topic;
}

double hdp_state_dynamic::compute_prior_probability_for_topic_sampling_online(
		int topic_id, bool decrease_table_counts) const {
	assert(topic_id <= m_num_topics);
	if (topic_id == m_num_topics) {
		return m_gamma;
	}
	int num_tables_from_prev_doc = 0;
	int actual_size = min((int) m_num_tables_by_zd_for_prev_doc_online.size(),
			m_num_topics);
	if (topic_id < actual_size) {
		num_tables_from_prev_doc =
				m_num_tables_by_zd_for_prev_doc_online[topic_id];
	}

	int num_tables_for_cur_doc = m_num_tables_by_zd_online[topic_id];
	int num_tables_by_z = m_num_tables_by_z[topic_id];
	if (decrease_table_counts) {
		--num_tables_for_cur_doc;
		--num_tables_by_z;
	}

	double prior_for_topic = num_tables_for_cur_doc + num_tables_from_prev_doc
			+ m_delta * num_tables_by_z;

	return prior_for_topic;
}

double hdp_state_dynamic::compute_normalisation_constant_for_topic_prior_online() const {
	return m_total_num_tables_by_d_online
			+ m_total_num_tables_by_d_for_prev_doc_online
			+ m_delta * m_total_num_tables + m_gamma;
}

void hdp_state_dynamic::sample_word_assignment_online(doc_state* d_state, int i,
		bool init, double_vec & q, double_vec & f) {
	if (init)
		doc_state_update_online(d_state, i, -1);

	if ((int) q.size() < d_state->m_num_tables + 1)
		q.resize(2 * d_state->m_num_tables + 1, 0.0);

	if ((int) f.size() < m_num_topics)
		f.resize(2 * m_num_topics + 1, 0.0);

	int k, t, w;
	w = d_state->m_words[i].m_word_index;
	double topic_new_prior =
			compute_prior_probability_for_topic_sampling_online(m_num_topics,
					false);
	double f_new = topic_new_prior / m_size_vocab;

	for (k = 0; k < m_num_topics; k++) {
		f[k] = (m_word_counts_by_zw[k][w] + m_eta)
				/ (m_word_counts_by_z[k] + m_size_vocab * m_eta);
		double topic_prior =
				compute_prior_probability_for_topic_sampling_online(k, false);
		f_new += topic_prior * f[k];
	}

	double prior_topic_normalisation =
			compute_normalisation_constant_for_topic_prior_online();
	f_new = f_new / prior_topic_normalisation;

	double total_q = 0.0, f_k = 0.0;
	for (t = 0; t < d_state->m_num_tables; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			k = d_state->m_table_to_topic[t];
			f_k = f[k];
		} else
			f_k = 0.0;

		total_q += d_state->m_word_counts_by_t[t] * f_k;
		q[t] = total_q;
	}
	total_q += m_alpha * f_new;
	q[d_state->m_num_tables] = total_q;

	double u = runiform() * total_q;
	for (t = 0; t < d_state->m_num_tables + 1; t++)
		if (u < q[t])
			break;

	d_state->m_words[i].m_table_assignment = t; // assign the new table

	if (t == d_state->m_num_tables) // this is a new table, we need get its k
			{
		f[m_num_topics] = 1.0 / m_size_vocab;
		int new_topic = sample_new_topic_for_new_table_online(f);
		doc_state_update_online(d_state, i, +1, new_topic);
	} else {
		doc_state_update_online(d_state, i, +1);
	}
}

// k is only provided when m_table_to_topic doesn't have that
void hdp_state_dynamic::doc_state_update(doc_state* d_state, int i, int update,
		int k) {
	int d, w, t;
	d = d_state->m_doc_id;
	w = d_state->m_words[i].m_word_index;
	//k = d_state->m_words[i].m_topic_assignment;
	t = d_state->m_words[i].m_table_assignment;
	if (k < 0)
		k = d_state->m_table_to_topic[t];
	assert(k >= 0);

	d_state->m_word_counts_by_t[t] += update;

	if (update == -1 && d_state->m_word_counts_by_t[t] == 0) /// this table becomes empty
			{
		m_total_num_tables--;
		m_total_num_tables_by_d[d]--;
		m_num_tables_by_zd[k][d]--;
		m_num_tables_by_z[k]--;
		d_state->m_table_to_topic[t] = -1;
		/// m_num_topics, no need to change at this moment
	}

	if (update == 1 && d_state->m_word_counts_by_t[t] == 1) /// a new table is created
			{
		if (t == d_state->m_num_tables)
			d_state->m_num_tables++; // create a new table
		d_state->m_table_to_topic[t] = k; // mapping the table

		if ((int) d_state->m_table_to_topic.size()
				< d_state->m_num_tables + 1) {
			d_state->m_table_to_topic.push_back(-1);
			d_state->m_word_counts_by_t.push_back(0);
		}
		if (k == m_num_topics) // used to k == m_num_topics
				{
			//assert(m_word_counts_by_z[k] == 1);
			if (k == m_num_topics)
				m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_zd.size() < m_num_topics + 1) {
				m_word_counts_by_z.push_back(0);
				m_num_tables_by_z.push_back(0);

				int* p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_word_counts_by_zd.push_back(p);

				p = new int[m_num_docs];
				memset(p, 0, sizeof(int) * m_num_docs);
				m_num_tables_by_zd.push_back(p);

				p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
		}
		m_num_tables_by_zd[k][d]++;          // adding the table to mixture k
		m_num_tables_by_z[k]++;
		m_total_num_tables_by_d[d]++;
		m_total_num_tables++;
	}

	m_word_counts_by_z[k] += update;
	m_word_counts_by_zw[k][w] += update;
	if (m_word_counts_by_zw[k][w] < 0) {
		cout << m_word_counts_by_zw[k][w];
	}
	m_word_counts_by_zd[k][d] += update;
}

// k is only provided when m_table_to_topic doesn't have that
void hdp_state_dynamic::doc_state_update_online(doc_state* d_state, int i,
		int update, int k) {
	int d, w, t;
	d = d_state->m_doc_id;
	w = d_state->m_words[i].m_word_index;
	//k = d_state->m_words[i].m_topic_assignment;
	t = d_state->m_words[i].m_table_assignment;
	if (k < 0)
		k = d_state->m_table_to_topic[t];
	assert(k >= 0);

	if (update == -1 && d_state->m_word_counts_by_t[t] == 0) /// this table becomes empty
			{
		m_total_num_tables_by_d_online--;
		m_total_num_tables--;
		m_num_tables_by_zd_online[k]--;
		m_num_tables_by_z[k]--;
		d_state->m_table_to_topic[t] = -1;

		/// m_num_topics, no need to change at this moment
	}

	d_state->m_word_counts_by_t[t] += update;

	if (update == 1 && d_state->m_word_counts_by_t[t] == 1) /// a new table is created
			{
		if (k == m_num_topics) // used to k == m_num_topics
				{
			//assert(m_word_counts_by_z[k] == 1);
			if (k == m_num_topics)
				m_num_topics++; // create a new topic
			if ((int) m_num_tables_by_zd.size() < m_num_topics + 1) {
				m_word_counts_by_z.push_back(0);
				m_num_tables_by_z.push_back(0);

				//int* p = new int [m_num_docs];
				//memset(p, 0, sizeof(int)*m_num_docs);
				//m_word_counts_by_zd.push_back(p);

				int * p = new int[m_size_vocab];
				memset(p, 0, sizeof(int) * m_size_vocab);
				m_word_counts_by_zw.push_back(p);
			}
			if ((int) m_word_counts_by_zd_online.size() < m_num_topics) {
				m_word_counts_by_zd_online.push_back(0);
				m_num_tables_by_zd_online.push_back(0);
			}
		}

		if (t == d_state->m_num_tables)
			d_state->m_num_tables++; // create a new table
		d_state->m_table_to_topic[t] = k; // mapping the table

		m_num_tables_by_zd_online[k]++;         // adding the table to mixture k
		m_total_num_tables_by_d_online++;
		m_num_tables_by_z[k]++;
		m_total_num_tables++;

		if ((int) d_state->m_table_to_topic.size()
				< d_state->m_num_tables + 1) {
			d_state->m_table_to_topic.push_back(-1);
			d_state->m_word_counts_by_t.push_back(0);
		}
	}

	m_word_counts_by_z[k] += update;
	m_word_counts_by_zw[k][w] += update;
	m_word_counts_by_zd_online[k] += update;
}

void hdp_state_dynamic::compact_doc_state(doc_state* d_state, int* k_to_new_k) {
	int num_tables_old = d_state->m_num_tables;
	int* t_to_new_t = new int[num_tables_old];

	int t, new_t, k, w;
	for (t = 0, new_t = 0; t < num_tables_old; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			t_to_new_t[t] = new_t;
			k = d_state->m_table_to_topic[t];
			d_state->m_table_to_topic[new_t] = k_to_new_k[k];
			swap_vec_element(d_state->m_word_counts_by_t, new_t, t);
			new_t++;
		} else
			d_state->m_table_to_topic[t] = -1;
	}
	d_state->m_num_tables = new_t;

	for (int i = 0; i < d_state->m_doc_length; i++) {
		t = d_state->m_words[i].m_table_assignment;
		new_t = t_to_new_t[t];
		d_state->m_words[i].m_table_assignment = new_t;
		w = d_state->m_words[i].m_word_index;
	}

	delete[] t_to_new_t;
}

void hdp_state_dynamic::compact_doc_state_online(doc_state* d_state,
		int* k_to_new_k) {
	int num_tables_old = d_state->m_num_tables;
	int* t_to_new_t = new int[num_tables_old];

	int t, new_t, k, w;
	for (t = 0, new_t = 0; t < num_tables_old; t++) {
		if (d_state->m_word_counts_by_t[t] > 0) {
			t_to_new_t[t] = new_t;
			k = d_state->m_table_to_topic[t];
			if (k >= m_num_topics_before_update) {
				d_state->m_table_to_topic[new_t] = k_to_new_k[k
						- m_num_topics_before_update];
			} else {
				d_state->m_table_to_topic[new_t] = k;
			}
			swap_vec_element(d_state->m_word_counts_by_t, new_t, t);
			new_t++;
		} else
			d_state->m_table_to_topic[t] = -1;
	}
	d_state->m_num_tables = new_t;

	for (int i = 0; i < d_state->m_doc_length; i++) {
		t = d_state->m_words[i].m_table_assignment;
		new_t = t_to_new_t[t];
		d_state->m_words[i].m_table_assignment = new_t;
		w = d_state->m_words[i].m_word_index;
	}

	delete[] t_to_new_t;
}

//compress the unused tables and components
void hdp_state_dynamic::compact_hdp_state() {
	int num_topics_old = m_num_topics;
	int* k_to_new_k = new int[num_topics_old];
	int k, new_k;
	for (k = 0, new_k = 0; k < num_topics_old; k++) {
		if (m_word_counts_by_z[k] > 0) {
			k_to_new_k[k] = new_k;
			swap_vec_element(m_word_counts_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_zd, new_k, k);
			swap_vec_element(m_word_counts_by_zd, new_k, k);
			swap_vec_element(m_word_counts_by_zw, new_k, k);
			new_k++;
		}
	}
	m_num_topics = new_k;

	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		d_state = m_doc_states[j];
		compact_doc_state(d_state, k_to_new_k);
	}

	delete[] k_to_new_k;
}

//compress the unused tables and components
void hdp_state_dynamic::compact_hdp_state_online() {
	int num_topics_old = m_num_topics;
	int* k_to_new_k = new int[m_num_topics - m_num_topics_before_update];
	int k, new_k;
	for (k = m_num_topics_before_update, new_k = m_num_topics_before_update;
			k < num_topics_old; k++) {
		if (m_word_counts_by_z[k] > 0) {
			k_to_new_k[k - m_num_topics_before_update] = new_k;
			swap_vec_element(m_word_counts_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_z, new_k, k);
			swap_vec_element(m_num_tables_by_zd_online, new_k, k);
			swap_vec_element(m_word_counts_by_zd_online, new_k, k);
			swap_vec_element(m_word_counts_by_zw, new_k, k);
			new_k++;
		}
	}
	m_num_topics = new_k;

	doc_state* d_state = NULL;
	for (int j = 0; j < m_num_docs; j++) {
		d_state = m_doc_states[j];
		compact_doc_state_online(d_state, k_to_new_k);
	}

	delete[] k_to_new_k;
}

void hdp_state_dynamic::save_state_to_file(char * name) {
	FILE * file = fopen(name, "wb");
	fwrite(&m_size_vocab, sizeof(int), 1, file);
	fwrite(&m_total_words, sizeof(int), 1, file);
	fwrite(&m_num_topics, sizeof(int), 1, file);
	fwrite(&m_total_num_tables_by_d[m_num_docs - 1], sizeof(int), 1, file);
	fwrite(&m_total_num_tables, sizeof(int), 1, file);

	fwrite(&m_eta, sizeof(double), 1, file);
	fwrite(&m_gamma, sizeof(double), 1, file);
	fwrite(&m_alpha, sizeof(double), 1, file);
	fwrite(&m_delta, sizeof(double), 1, file);

	for (int k = 0; k < m_num_topics; k++) {
		fwrite(&(m_num_tables_by_zd[k][m_num_docs - 1]), sizeof(int), 1, file);
		fwrite(&(m_num_tables_by_z[k]), sizeof(int), 1, file);
		fwrite(&(m_word_counts_by_z[k]), sizeof(int), 1, file);
		fwrite(m_word_counts_by_zw[k], sizeof(int), m_size_vocab, file);
	}
	fclose(file);
}

void hdp_state_dynamic::save_state_ex(char * name, int sample_num) {

	if (sample_num == -1) {
		save_state_to_file(name);
	} else {
		char * file_name = new char[strlen(name) + 5 * sizeof(char)];
		strcpy(file_name, name);
		char* postfix = new char[5];
		postfix[0] = '_';
		postfix[1] = '0' + sample_num / 100;
		postfix[2] = '0' + (sample_num / 10) % 10;
		postfix[3] = '0' + sample_num % 10;
		postfix[4] = '\0';

		file_name = strcat(file_name, postfix);
		save_state_to_file(file_name);
		delete[] file_name;
		delete[] postfix;
	}
}

void hdp_state_dynamic::load_state_from_file(char * name) {
	free_state();
	FILE * file = fopen(name, "rb");
	fread(&m_size_vocab, sizeof(int), 1, file);
	fread(&m_total_words, sizeof(int), 1, file);
	fread(&m_num_topics, sizeof(int), 1, file);
	fread(&m_total_num_tables_by_d_for_prev_doc_online, sizeof(int), 1, file);
	fread(&m_total_num_tables, sizeof(int), 1, file);

	fread(&m_eta, sizeof(double), 1, file);
	fread(&m_gamma, sizeof(double), 1, file);
	fread(&m_alpha, sizeof(double), 1, file);
	fread(&m_delta, sizeof(double), 1, file);

	m_num_tables_by_zd_for_prev_doc_online.resize(m_num_topics);
	m_word_counts_by_z.resize(m_num_topics);
	m_num_tables_by_z.resize(m_num_topics);
	m_word_counts_by_zw.resize(m_num_topics);
	for (int k = 0; k < m_num_topics; k++) {
		fread(&(m_num_tables_by_zd_for_prev_doc_online[k]), sizeof(int), 1,
				file);
		fread(&(m_num_tables_by_z[k]), sizeof(int), 1, file);
		fread(&(m_word_counts_by_z[k]), sizeof(int), 1, file);

		m_word_counts_by_zw[k] = new int[m_size_vocab];
		fread(m_word_counts_by_zw[k], sizeof(int), m_size_vocab, file);
	}
	fclose(file);

}

void hdp_state_dynamic::load_state_ex(char * name, int sample_num) {
	if (sample_num == -1) {
		load_state_from_file(name);
	} else {
		char * file_name = new char[strlen(name) + 5 * sizeof(char)];
		strcpy(file_name, name);
		char* postfix = new char[5];
		postfix[0] = '_';
		postfix[1] = '0' + sample_num / 100;
		postfix[2] = '0' + (sample_num / 10) % 10;
		postfix[3] = '0' + sample_num % 10;
		postfix[4] = '\0';

		file_name = strcat(file_name, postfix);
		load_state_from_file(file_name);
		delete[] file_name;
		delete[] postfix;
	}
}

void hdp_state_dynamic::save_online_update_for_next_iteration() {
	m_num_tables_by_zd_for_prev_doc_online = m_num_tables_by_zd_online;
	m_total_num_tables_by_d_for_prev_doc_online =
			m_total_num_tables_by_d_online;
}

void hdp_state_dynamic::update_state_after_online_update() {
	save_online_update_for_next_iteration();
	free_doc_info();
	free_current_doc_counts_update();
}

void hdp_state_dynamic::free_current_doc_counts_update() {
	m_num_tables_by_zd_online.clear();
	m_word_counts_by_zd_online.clear();
	m_total_num_tables_by_d_online = 0;
}

void hdp_state_dynamic::save_last_document_info_for_online_update() {
	m_num_tables_by_zd_for_prev_doc_online.clear();
	m_num_tables_by_zd_for_prev_doc_online.assign(m_num_topics, 0);

	for (int k = 0; k < m_num_topics; ++k) {
		m_num_tables_by_zd_for_prev_doc_online[k] =
				m_num_tables_by_zd[k][m_num_docs - 1];
	}

	m_total_num_tables_by_d_for_prev_doc_online =
			m_total_num_tables_by_d[m_num_docs - 1];
}

