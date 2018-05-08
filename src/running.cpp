#include "running.h"
#include "corpus.h"
#include "hdp.h"
#include "utils.h"
#include "likelihood_computer.h"
#include <string>
#include <fstream>
using namespace std;



void run_train_batch_whole_process(int m_samplers, int vocab_size, double gamma_a, double gamma_b, double alpha_a, 
								   double alpha_b, double eta, double gamma, double alpha,
                                   int max_iter,
								   const char * directory, const char * data_path, char * model_path,
								   bool save_topic_distribution)
{
	// read data
    corpus * c = new corpus();
    c->read_data(data_path);

    // read hyperparameters

    hdp_hyperparameter * hdp_hyperparam = new hdp_hyperparameter();
    hdp_hyperparam->setup_parameters(gamma_a, gamma_b, alpha_a, alpha_b, max_iter);

    hdp ** hdp_instance = new hdp * [m_samplers];
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s] = new hdp();
    }

	for (int s = 0; s < m_samplers; ++s) {
		printf("batch run for sampler = %05d\n", s);
		hdp_instance[s]->setup_state(c, eta, gamma, alpha, *hdp_hyperparam, vocab_size);

        hdp_instance[s]->run();

		if (model_path != NULL) {
			hdp_instance[s]->save(model_path, s);
		}
		if (save_topic_distribution) {
			string file_name = directory;
            file_name += "/topic_distributions";
            file_name = add_postfix(file_name, s);
            print_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd,
                         hdp_instance[s]->m_state->m_num_docs,
                         file_name);
		}
	}

    // free resources
    if (hdp_instance != NULL) {
        for (int s = 0; s < m_samplers; ++s) {
            delete hdp_instance[s];
        }
    }
    delete [] hdp_instance;
    delete hdp_hyperparam;
	delete c;
}

void run_online_process(int m_samplers, int vocab_size, double gamma_a,
                        double gamma_b, double alpha_a, double alpha_b,
						double eta, double gamma, double alpha, int max_iter,
                        const char * directory, const char * data_path_batch,
						const char * data_path_online, char * model_path,
						bool hm_likelihood_computing,
						bool save_topic_distribution,
                        int _m_online_documents)
{
	hdp ** hdp_instance = new hdp * [m_samplers];
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s] = new hdp();
    }

	
	// read hyperparameters

	hdp_hyperparameter * hdp_hyperparam = new hdp_hyperparameter();
	hdp_hyperparam->setup_parameters(gamma_a, gamma_b,
										alpha_a, alpha_b,
										max_iter);

	if (model_path == NULL) {
		// read data
		corpus * c = new corpus();
		c->read_data(data_path_batch);

        for (int s = 0; s < m_samplers; ++s) {
            printf("batch run for sampler = %05d, ", s);
            hdp_instance[s]->setup_state(c, eta, gamma, alpha,
                                         *hdp_hyperparam, vocab_size);

            hdp_instance[s]->run();

			if (save_topic_distribution) {
				string file_name = directory;
				file_name += "/topic_distributions";
				file_name = add_postfix(file_name, s);
				print_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd,
							 hdp_instance[s]->m_state->m_num_docs,
							 file_name);
			}
        }

		delete c;
	} else {
        for (int s = 0; s < m_samplers; ++s) {
            hdp_instance[s]->load(model_path, s);
            hdp_instance[s]->setup_state(*hdp_hyperparam);
        }
	}
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s]->free_online_info();
		
		if (save_topic_distribution) {
			string file_name = directory;
			file_name += "/online_topic_distributions";
			file_name = add_postfix(file_name, s);
			ofstream out_fst;
			out_fst.open(file_name);
			out_fst.close();
		}

    }

	FILE * fileptr;

	fileptr = fopen(data_path_online, "r");

	bool eof_online = false;
    likelihood_computer l_computer;
	counts ** counts_before_update = new counts * [m_samplers];
	for (int s  = 0; s < m_samplers; ++s) {
		counts_before_update[s] = new counts;
	}
	
	int m_num_online_docs = 0;
    
	document * doc = new document();
    doc->read_data(fileptr, eof_online, m_num_online_docs);
    
    vector<double> hm_likelihood_vec;
    hm_likelihood_vec.clear();
    
    if (hm_likelihood_computing) {
        if (_m_online_documents > 0) {
            hm_likelihood_vec.assign(_m_online_documents, 0);
        }
    }
	
	while (!eof_online) {
		++m_num_online_docs;

        // compute online posterior estimates
		if (doc->length > 0) {
			for (int s = 0; s < m_samplers; ++s) {
				printf("online doc = %05d, sampler = %05d, ", m_num_online_docs, s);

				if (hm_likelihood_computing) {
					counts_before_update[s]->set_counts(hdp_instance[s]->m_state->m_size_vocab,
														hdp_instance[s]->m_state->m_total_num_tables,
														hdp_instance[s]->m_state->m_num_tables_by_z,
														hdp_instance[s]->m_state->m_word_counts_by_z,
														hdp_instance[s]->m_state->m_word_counts_by_zw);
				}

				hdp_instance[s]->setup_doc_info_update(doc);

                hdp_instance[s]->run_online();

				if (save_topic_distribution) {
					string file_name = directory;
					file_name += "/online_topic_distributions";
					file_name = add_postfix(file_name, s);
					append_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd_online,
									file_name);
				}
			}
		}

        // compute predictive likelihood via harmonic mean
		if (hm_likelihood_computing) {
			double hm_likelihood = 0;
			if (doc->length > 0) {
				hm_likelihood = l_computer.compute_harmonic_mean_predictive_likelihood(hdp_instance, counts_before_update, m_samplers);
			}

			printf("online doc = %05d, harmonic mean likelihood = %.5f\n", m_num_online_docs, hm_likelihood);
            
            if (hm_likelihood_vec.size() >= m_num_online_docs) {
                hm_likelihood_vec[m_num_online_docs - 1] = hm_likelihood / doc->length;
            } else {
                hm_likelihood_vec.push_back(hm_likelihood / doc->length);
            }
		}

		// free doc info
		if (doc->length > 0) {
			for (int s = 0; s < m_samplers; ++s) {
				hdp_instance[s]->free_online_info();
			}
		}
            
        doc->read_data(fileptr, eof_online, m_num_online_docs);
	}
    delete doc;
	fclose(fileptr); // close the file

    if (hm_likelihood_computing) {
        string file_name = directory;
        file_name += "/predictive_hm_likelihood";
        print_vector(hm_likelihood_vec, file_name);
    }
    
    // free resources
    if (hdp_instance != NULL) {
        for (int s = 0; s < m_samplers; ++s) {
            delete hdp_instance[s];
        }
		delete [] hdp_instance;
    }
    delete hdp_hyperparam;

	if (counts_before_update != NULL) {
		for (int s  = 0; s < m_samplers; ++s) {
			delete counts_before_update[s];
		}
		delete [] counts_before_update;
	}
}


void run_train_dynamic_batch_whole_process(int m_samplers, int vocab_size, double gamma_a,
                                           double gamma_b, double alpha_a, double alpha_b,
                                           double eta, double gamma, double alpha,
                                           double delta, int max_iter,
                                           const char * directory, const char * data_path,
                                           char * model_path,
                                           bool save_topic_distributions)
{
    // read data
    corpus * c = new corpus();
    c->read_data(data_path);
    
    // read hyperparameters
    
    hdp_hyperparameter * hdp_hyperparam = new hdp_hyperparameter();
    hdp_hyperparam->setup_parameters(gamma_a, gamma_b,
                                     alpha_a, alpha_b, max_iter);
    
    hdp_dynamic ** hdp_instance = new hdp_dynamic * [m_samplers];
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s] = new hdp_dynamic();
    }
    
    for (int s = 0; s < m_samplers; ++s) {
        printf("batch run for sampler = %05d\n", s);
        hdp_instance[s]->setup_state(c, eta, gamma, alpha, delta, *hdp_hyperparam, vocab_size);

        hdp_instance[s]->run();
        
        if (!(model_path == NULL)) {
            hdp_instance[s]->save(model_path, s);
        }

        if (save_topic_distributions) {
            string file_name = directory;
            file_name += "/topic_distributions";
            file_name = add_postfix(file_name, s);
            print_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd,
                         hdp_instance[s]->m_state->m_num_docs,
                         file_name);
        }
    }
    
    // free resources
    if (hdp_instance != NULL) {
        for (int s = 0; s < m_samplers; ++s) {
            delete hdp_instance[s];
        }
    }
    delete [] hdp_instance;
    delete hdp_hyperparam;
    delete c;
}

void run_online_dynamic_process(int m_samplers, int vocab_size,
                                double gamma_a, double gamma_b, double alpha_a,
                                double alpha_b, double eta, double gamma, double alpha,
                                double delta, int max_iter, const char * directory,
                                const char * data_path_batch,
                                const char * data_path_online, char * model_path,
                                bool hm_likelihood_computing,
								bool save_topic_distribution,
                                int _m_online_documents)
{
    hdp_dynamic ** hdp_instance = new hdp_dynamic * [m_samplers];
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s] = new hdp_dynamic();
    }
    
    
    // read hyperparameters
    
    hdp_hyperparameter * hdp_hyperparam = new hdp_hyperparameter();
    hdp_hyperparam->setup_parameters(gamma_a, gamma_b,
                                     alpha_a, alpha_b,
                                     max_iter);
    
    if (model_path == NULL) {
        // read data
        corpus * c = new corpus();
        c->read_data(data_path_batch);
        
        for (int s = 0; s < m_samplers; ++s) {
            printf("batch run for sampler = %05d, ", s);
            hdp_instance[s]->setup_state(c, eta, gamma, alpha, delta, *hdp_hyperparam, vocab_size);

            hdp_instance[s]->run();
            
            hdp_instance[s]->save_last_document_info_for_online_update();

			if (save_topic_distribution) {
				string file_name = directory;
				file_name += "/topic_distributions";
				file_name = add_postfix(file_name, s);
				print_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd,
							 hdp_instance[s]->m_state->m_num_docs,
							 file_name);
			}
        }
        
        delete c;
    } else {
        for (int s = 0; s < m_samplers; ++s) {
            hdp_instance[s]->load(model_path, s);
            hdp_instance[s]->setup_state(*hdp_hyperparam);
        }
    }
    for (int s = 0; s < m_samplers; ++s) {
        hdp_instance[s]->free_doc_info();
		if (save_topic_distribution) {
			string file_name = directory;
			file_name += "/online_topic_distributions";
			file_name = add_postfix(file_name, s);
			ofstream out_fst;
			out_fst.open(file_name);
			out_fst.close();
		}
    }
    
    FILE * fileptr;
    
    fileptr = fopen(data_path_online, "r");
    
    bool eof_online = false;

    likelihood_computer l_computer;
    counts ** counts_before_update = new counts * [m_samplers];
    for (int s  = 0; s < m_samplers; ++s) {
        counts_before_update[s] = new counts;
    }
    
    int m_num_online_docs = 0;
    
    document * doc = new document();
    doc->read_data(fileptr, eof_online, m_num_online_docs);
    
    vector<double> hm_likelihood_vec;
    hm_likelihood_vec.clear();
    
    if (hm_likelihood_computing) {
        if (_m_online_documents > 0) {
            hm_likelihood_vec.assign(_m_online_documents, 0);
        }
    }
    
    while (!eof_online) {
        ++m_num_online_docs;
        
        // compute online posterior estimates
		if (doc->length > 0) {
			for (int s = 0; s < m_samplers; ++s) {
				printf("online doc = %05d, sampler = %05d, ", m_num_online_docs - 1, s);
            
				if (hm_likelihood_computing) {
					counts_before_update[s]->set_counts(hdp_instance[s]->m_state->m_size_vocab,
														hdp_instance[s]->m_state->m_total_num_tables,
														hdp_instance[s]->m_state->m_num_tables_by_z,
														hdp_instance[s]->m_state->m_word_counts_by_z,
														hdp_instance[s]->m_state->m_word_counts_by_zw);
				}
            
				hdp_instance[s]->setup_doc_info_update(doc);

                hdp_instance[s]->run_online();


				if (save_topic_distribution) {
					string file_name = directory;
					file_name += "/online_topic_distributions";
					file_name = add_postfix(file_name, s);
					append_matrix(hdp_instance[s]->m_state->m_word_counts_by_zd_online,
								  file_name);
				}
			}
		}
        
        // compute predictive likelihood via harmonic mean
        if (hm_likelihood_computing) {
            double hm_likelihood = 0;
			if (doc->length > 0) {
				hm_likelihood = 
					l_computer.compute_harmonic_mean_predictive_likelihood(hdp_instance, counts_before_update, m_samplers);
			}
            
			printf("online doc = %05d, harmonic mean likelihood = %.5f\n", m_num_online_docs, hm_likelihood);
            
            if (hm_likelihood_vec.size() >= m_num_online_docs) {
                hm_likelihood_vec[m_num_online_docs - 1] = hm_likelihood / doc->total;
            } else {
                hm_likelihood_vec.push_back(hm_likelihood / doc->total);
            }
        }
        
        
        // free doc info
		if (doc->length > 0) {
			for (int s = 0; s < m_samplers; ++s) {
				hdp_instance[s]->update_state_after_online_update();
			}
		}
        
        doc->read_data(fileptr, eof_online, m_num_online_docs);
    }
    delete doc;
    fclose(fileptr); // close the file

    if (hm_likelihood_computing) {
        string file_name = directory;
        file_name += "/predictive_hm_likelihood";
        print_vector(hm_likelihood_vec, file_name);
    }
    
    // free resources
    if (hdp_instance != NULL) {
        for (int s = 0; s < m_samplers; ++s) {
            delete hdp_instance[s];
        }
        delete [] hdp_instance;
    }
    delete hdp_hyperparam;
    
    if (counts_before_update != NULL) {
        for (int s  = 0; s < m_samplers; ++s) {
            delete counts_before_update[s];
        }
        delete [] counts_before_update;
    }
}

