#ifndef STATE_H
#define STATE_H

#include "corpus.h"
#include <map>
#include <vector>
using namespace std;

class hdp_hyperparameter
{
    /// hyperparameters
public:
    double m_gamma_a;
    double m_gamma_b;
    double m_alpha_a;
    double m_alpha_b;
    int    m_max_iter;

public:
    void setup_parameters(double _gamma_a, double _gamma_b,
                        double _alpha_a, double _alpha_b,
                        int _max_iter)
    {
        m_gamma_a   = _gamma_a;
        m_gamma_b   = _gamma_b;
        m_alpha_a   = _alpha_a;
        m_alpha_b   = _alpha_b;
        m_max_iter  = _max_iter;
    }
    void copy_parameters(const hdp_hyperparameter _input_hyperparam);
};

typedef vector<int> int_vec; // define the vector of int
typedef vector<double> double_vec; // define the vector of double

/// word info structure used in the main class
struct word_info
{
public:
    int m_word_index;
    int m_table_assignment;
    //int m_topic_assignment; // this is extra information
};

class doc_state
{
public:
    int m_doc_id; // document id
    int m_doc_length;  // document length
    int m_num_tables;  // number of tables in this document
    word_info * m_words;

    int_vec m_table_to_topic; // for a doc, translate its table index to topic index
    int_vec m_word_counts_by_t; // word counts for each table
    
    //vector < vector<int> > m_words_by_zi; // stores the word idx indexed by z then i
public:
    doc_state();
    ~doc_state();
public:
    void setup_state_from_doc(const document * doc);
    void free_doc_state();
};

class counts
{
public:
    int m_num_words;
	int m_total_num_tables;
    int_vec m_num_tables_by_z;
    int_vec m_word_counts_by_z;
    vector <int*> m_word_counts_by_zw;
    
    counts();
    counts(int _m_num_words, int _m_total_num_tables, const int_vec& _m_num_tables_by_z, const int_vec& _m_word_counts_by_z, const vector<int*>& _m_word_counts_by_zw);
    ~counts();

	void set_counts(int _m_num_words, int _m_total_num_tables, const int_vec& _m_num_tables_by_z,
                    const int_vec& _m_word_counts_by_z,
                    const vector<int*>& _m_word_counts_by_zw);
};

class word_counts
{
public:
    int_vec m_word_counts_by_z;
    vector<int*> m_word_counts_by_zw;
    
    word_counts();
    ~word_counts();
};

class hdp_state
{
public:

/// doc information, fix value
    int m_size_vocab;
    int m_total_words;
    int m_num_docs;

/// document states
    doc_state** m_doc_states;

/// number of topics
    int m_num_topics;


/// total number of tables for all topics
    int m_total_num_tables;

/// by_z, by topic
/// by_d, by document, for each topic
/// by_w, by word, for each topic
/// by_t, by table for each document
    int_vec   m_num_tables_by_z; // how many tables each topic has
    int_vec   m_word_counts_by_z;   // word counts for each topic
    vector <int*> m_word_counts_by_zd; // word counts for [each topic, each doc]
    vector <int*> m_word_counts_by_zw; // word counts for [each topic, each word]

/// for online updates (local changes for the current document)
//	int_vec m_num_tables_by_z_update;
//	int_vec m_word_counts_by_z_update;
	int_vec m_word_counts_by_zd_online;
//	vector <int*> m_word_counts_by_zw_update;
	int m_num_topics_before_update;	// number of topics from the batch plus possible new topics
    
///

/// topic Dirichlet parameter
    double m_eta;

/// including concentration parameters
    double m_gamma;
    double m_alpha;
public:
    hdp_state();
    virtual ~hdp_state();
public:
	void   set_vocab_size(int _vocab_size);
    void   setup_state_from_corpus(const corpus* c);
	void   setup_doc_info_from_document(const document * doc);
    void   allocate_initial_space();
	void   allocate_space_for_online_updates();
    void   free_state();
	void   free_doc_info();
    void   free_counts_update();
    void   free_online_info();

    void   iterate_gibbs_state(bool remove, bool permute,
                               hdp_hyperparameter* hdp_hyperparam,
                               bool table_sampling=false);
	void   iterate_gibbs_state_online(bool remove, hdp_hyperparameter* hdp_hyperparam,
									  bool table_sampling=false);

    void   sample_tables(doc_state* d_state, double_vec & q, double_vec & f);
	void   sample_tables_online(doc_state* d_state, double_vec & q, double_vec & f);

    void   sample_table_assignment(doc_state* d_state, int t, int* words, double_vec & q, double_vec & f);
	void   sample_table_assignment_online(doc_state* d_state, int t, int* words, double_vec & q, double_vec & f);

    void   sample_word_assignment(doc_state* d_state, int i, bool remove, double_vec & q, double_vec & f);
	void   sample_word_assignment_online(doc_state* d_state, int i, bool remove, double_vec & q, double_vec & f);

    void   doc_state_update(doc_state* d_state, int i, int update, int k=-1);
	void   doc_state_update_online(doc_state* d_state, int i, int update, int k=-1);

    void   compact_doc_state(doc_state* d_state, int* k_to_new_k);
	void   compact_doc_state_online(doc_state* d_state, int* k_to_new_k);

    void   compact_hdp_state();
	void   compact_hdp_state_online();

	void   save_state_to_file(char * name);
    void   save_state(char * name);
    void   save_state_ex(char * name, int sampler_num = -1);
	void   load_state_from_file(char * name);
    void   load_state_ex(char * name, int sampler_num = -1);
};

class hdp_state_dynamic
{
public:

/// doc information, fix value
    int m_size_vocab;
    int m_total_words;
    int m_num_docs;

/// document states
    doc_state** m_doc_states;

/// number of topics
    int m_num_topics;

/// total number of tables for all topics (during the iteration of batch this is the total number of tables for the documents seen so far)
    int m_total_num_tables;


/// by_z, by topic
/// by_d, by document, for each topic
/// by_w, by word, for each topic
/// by_t, by table for each document
 
	int_vec   m_word_counts_by_z;   // word counts for each topic
	int_vec   m_num_tables_by_z;	// how many tables each topic has (during the iterations of batch this is the number of tables for the documents seen so far)
	vector<int*> m_num_tables_by_zd; // how many tables [each topic, (in) each doc] has
	int_vec m_total_num_tables_by_d; // total number of tables for each doc
    vector <int*> m_word_counts_by_zd; // word counts for [each topic, each doc]
    vector <int*> m_word_counts_by_zw; // word counts for [each topic, each word]

/// for online updates (local changes for the current document)
	int_vec m_word_counts_by_zd_online;
	int_vec m_num_tables_by_zd_online;
	int_vec m_num_tables_by_zd_for_prev_doc_online;
	int m_total_num_tables_by_d_online;
	int m_total_num_tables_by_d_for_prev_doc_online;
	int m_num_topics_before_update;	// number of topics from the batch
    
///

/// topic Dirichlet parameter
    double m_eta;

/// including concentration parameters
    double m_gamma;
    double m_alpha;
	double m_delta;
public:
    hdp_state_dynamic();
    virtual ~hdp_state_dynamic();
public:
	void   set_vocab_size(int _vocab_size);
    void   setup_state_from_corpus(const corpus* c);
	void   setup_doc_info_from_document(const document * doc);
    void   allocate_initial_space();
	void   allocate_space_for_online_updates();

    void   free_state();
	void   free_doc_info();
    void   free_counts_update();
    void   free_online_info();
    void   free_current_doc_counts_update();

    void   iterate_gibbs_state(bool remove,
                               hdp_hyperparameter* hdp_hyperparam,
                               bool table_sampling=false);
	void   iterate_gibbs_state_online(bool remove, hdp_hyperparameter* hdp_hyperparam,
									  bool table_sampling=false);

	void   reset_counts_up_to_date(int doc_id = 0);

    void   sample_tables(doc_state* d_state, bool not_init, double_vec & q, double_vec & f);
	void   sample_tables_online(doc_state* d_state, double_vec & q, double_vec & f);

    void   sample_table_assignment(doc_state* d_state, int t, bool not_init,
                                   int* words, double_vec & q, double_vec & f);
	void   sample_table_assignment_online(doc_state* d_state, int t, int* words, double_vec & q, double_vec & f);

    void   sample_word_assignment(doc_state* d_state, int i, bool remove, double_vec & q, double_vec & f);
	void   sample_word_assignment_online(doc_state* d_state, int i, bool remove, double_vec & q, double_vec & f);

	int    sample_new_topic_for_new_table_online(const vector<double>& likelihood_per_topic_assignment) const;

	double compute_prior_probability_for_topic_sampling(int doc_id, int topic_id,
                                                        bool decrease_table_counts = false) const;
	double compute_prior_probability_for_topic_sampling_online(int topic_id,
                                                               bool decrease_table_counts = false) const;

	double compute_normalisation_constant_for_topic_prior_online() const;

	void   doc_state_update(doc_state* d_state, int i, int update, int k=-1);
	void   doc_state_update_online(doc_state* d_state, int i, int update, int k=-1);

    void   compact_doc_state(doc_state* d_state, int* k_to_new_k);
	void   compact_doc_state_online(doc_state* d_state, int* k_to_new_k);

    void   compact_hdp_state();
	void   compact_hdp_state_online();

	double compute_next_doc_table_assignment_likelihood(int doc_id, int old_topic_id = -1, int new_topic_id = -1, bool new_table = false) const;
    double compute_next_table_assignment_likelihood(int doc_id, int old_topic_id = -1,
                                                    int new_topic_id = -1, bool new_table = false) const;
    double compute_next_table_assignment_likelihood_for_non_born_topics(int doc_id,
                                                                        const vector<int>& updating_num_tables_by_z,
                                                                        int old_topic_id = -1,
                                                                        int new_topic_id = -1,
                                                                        bool update_counts_for_current_doc = false) const;

	void   save_state_to_file(char * name);
    void   save_state_ex(char * name, int sampler_num = -1);
	void   load_state_from_file(char * name);
    void   load_state_ex(char * name, int sampler_num = -1);
    
    void save_online_update_for_next_iteration();
    void update_state_after_online_update();
    void save_last_document_info_for_online_update();
};

#endif // STATE_H
