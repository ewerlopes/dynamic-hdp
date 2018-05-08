#ifndef RUNNING_H
#define RUNNING_H


void run_train_batch_whole_process(int m_samplers, int vocab_size, double gamma_a, double gamma_b, 
								   double alpha_a, double alpha_b, double eta, double gamma, double alpha,
								   int max_iter,
								   const char * directory,
								   const char * data_path, char * model_path,
								   bool save_topic_distribution);

void run_online_process(int m_samplers, int vocab_size, double gamma_a, double gamma_b,
                        double alpha_a, double alpha_b,
						double eta, double gamma, double alpha, int max_iter,
						const char * directory, const char * data_path_batch,
						const char * data_path_online, char * model_path,
						bool hm_likelihood, bool save_topic_distribution, int _m_online_documents = -1);

void run_train_dynamic_batch_whole_process(int m_samplers, int vocab_size, double gamma_a,
                                           double gamma_b, double alpha_a, double alpha_b,
                                           double eta, double gamma, double alpha,
										   double delta, int max_iter, const char * directory,
                                           const char * data_path, char * model_path,
										   bool save_topic_distributions);

void run_online_dynamic_process(int m_samplers, int vocab_size,
                                double gamma_a, double gamma_b, double alpha_a,
                                double alpha_b, double eta, double gamma, double alpha,
								double delta, int max_iter, const char * directory,
                                const char * data_path_batch, const char * data_path_online,
                                char * model_path, bool hm_likelihood,
								bool save_topic_distribution,
                                int _m_online_documents = -1);

#endif // RUNNING_H