#ifndef HDP_H
#define HDP_H

#include "state.h"


/// implement the Chinese restaurant franchies algorithm
class hdp
{
public:
/// fixed parameters
    hdp_hyperparameter * m_hdp_param;

/// sampling state
    hdp_state * m_state;

public:
    hdp();
    virtual ~hdp();
public:
    void run();
	void run_online();

	void setup_state(const corpus * c,
                     double _eta, double _gamma, double _alpha,
                     const hdp_hyperparameter _hdp_param, int vocab_size = 0);

	void setup_state(const hdp_hyperparameter _hdp_param);
	void setup_doc_info_update(const document * doc);

	int get_size_vocab() const;
	double get_eta() const; 
	int get_topics_number() const;
    int get_word_topic_counts(int k) const;
    int get_word_topic_counts(int k, int w) const;

	void load(char * model_path, int sample_num = -1);
	void save(char * model_path, int sample_num = -1);

	void free_online_info();
    void free_hdp();
};

class hdp_dynamic
{
public:
/// fixed parameters
    hdp_hyperparameter * m_hdp_param;

/// sampling state
    hdp_state_dynamic * m_state;

public:
    hdp_dynamic();
    virtual ~hdp_dynamic();
public:
    void run();

	void run_online();

	void setup_state(const corpus *c, double _eta, double _gamma, double _alpha, double _delta,
                     const hdp_hyperparameter _hdp_param, int vocab_size);

	void setup_state(const hdp_hyperparameter _hdp_param);
	void setup_doc_info_update(const document * doc);

	int get_size_vocab() const;
	double get_eta() const; 
	int get_topics_number() const;
    int get_word_topic_counts(int k) const;
    int get_word_topic_counts(int k, int w) const;

	void load(char * model_path, int sample_num = -1);
	void save(char * model_path, int sample_num = -1);


	void free_doc_info();
    void free_hdp();

	void update_state_after_online_update();
	void save_last_document_info_for_online_update();
};


#endif // HDP_H
