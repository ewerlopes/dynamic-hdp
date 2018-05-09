#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "utils.h"
#include "hdp.h"
#include "hdp_dynamic.h"
#include "running.h"
#define VERBOSE true

gsl_rng * RANDOM_NUMBER;

void print_usage_and_exit()
{
    printf("\nC++ implementation of Gibbs sampling for dynamic hierarchical Dirichlet process.\n");
    printf("Authors: {o.isupova, dkuzin1}@sheffield.ac.uk, Department of Automatic Control and Systems Engineering, University of Sheffield.\n");
    printf("usage:\n");
    printf("      hdp [options]\n");
    printf("      general parameters:\n");
    printf("      --algorithm:      train/online/train_dynamic/online_dynamic, not optional\n");
    printf("      --data:           train data file, in lda-c format, not optional\n");
    printf("      --online_data:    testing data file, in lda-c format, not optional for online inference\n");
    printf("      --directory:      save directory\n");
    printf("      --max_iter:       the max number of iterations, default 1000\n");
    printf("      --random_seed:    the random seed, default from the current time\n");

    printf("      --sampler_number: the number of posterior samples for computing likelihood, default 1\n");
    printf("      --vocab_size:     vocabulary size\n");
    printf("      --online_data_size: testing data size\n");

    printf("\n      training parameters:\n");
    printf("      --gamma_a:        shape for 1st-level concentration parameter, default 1.0\n");
    printf("      --gamma_b:        scale for 1st-level concentration parameter, default 1.0\n");
    printf("      --alpha_a:        shape for 2nd-level concentration parameter, default 1.0\n");
    printf("      --alpha_b:        scale for 2nd-level concentration parameter, default 1.0\n");
    printf("      --delta:          weight for the global topic usage factor, default 1.0\n");
    printf("      --eta:            topic Dirichlet parameter, default 0.5\n");

    printf("      --model_path:    path for saved model\n");
    printf("      --hm_likelihood_online: compute likelihood for testing data yes or no, default no\n");

    printf("\nexamples:\n");
    printf("      ./hdp --algorithm train --data data --directory train_dir --model_path model_path --vocab_size vocabulary_size --sampler_number number_of_samples\n");
    printf("      ./hdp --algorithm online --data data --online_data test_data --model_path saved_model --directory test_dir --vocab_size vocabulaty_size --sampler_number number_of_samples --hm_likelihood_online yes\n");
    printf("      ./hdp --algorithm train_dynamic --data data --directory train_dir --model_path model_path --vocab_size vocabulary_size --sampler_number number_of_samples\n");
    printf("      ./hdp --algorithm online_dynamic --data data --online_data test_data --model_path saved_model --directory test_dir --vocab_size vocabulaty_size --sampler_number number_of_samples --hm_likelihood_online yes\n");
    printf("\n");

	char z;
	cout << ">>";
	cin >> z;
    exit(0);
}

int main(int argc, char** argv)
{
    if (argc < 2 || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help") ||
            !strcmp(argv[1], "-h") || !strcmp(argv[1], "--usage"))
    {
        print_usage_and_exit();
    }

    double gamma_a = 1.0;
	double gamma_b = 1.0;
	double alpha_a = 1.0;
	double alpha_b = 1.0;
	double eta = 0.5;
	double delta = 1;
    double gamma = -1;
    double alpha = -1;
    int max_iter = 1000;
	bool harm_mean_likelihood_online = false;
    bool save_topic_distribution = true;
	int m_vocab_size = 0;
    int m_samplers = 1;
    int m_online_documents = -1;

    time_t t;
    time(&t);
    long seed = (long) t;

    char* directory = NULL;
    char* algorithm = NULL;
    char* data_path = NULL;
	char* data_path_online = NULL;
    char* model_path = NULL;

    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "--algorithm"))        algorithm = argv[++i];
        else if (!strcmp(argv[i], "--data"))        data_path = argv[++i];
		else if (!strcmp(argv[i], "--online_data")) data_path_online = argv[++i];
        else if (!strcmp(argv[i], "--max_iter"))    max_iter = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--directory"))   directory = argv[++i];
        else if (!strcmp(argv[i], "--random_seed")) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--sampler_number")) m_samplers = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--vocab_size"))  m_vocab_size = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--online_data_size"))  m_online_documents = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--gamma_a"))     gamma_a = atof(argv[++i]);
        else if (!strcmp(argv[i], "--gamma_b"))     gamma_b = atof(argv[++i]);
        else if (!strcmp(argv[i], "--alpha_a"))     alpha_a = atof(argv[++i]);
        else if (!strcmp(argv[i], "--alpha_b"))     alpha_b = atof(argv[++i]);
        else if (!strcmp(argv[i], "--eta"))         eta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--gamma"))       gamma = atof(argv[++i]);
        else if (!strcmp(argv[i], "--alpha"))       alpha = atof(argv[++i]);
        else if (!strcmp(argv[i], "--delta"))       delta = atof(argv[++i]);
        else if (!strcmp(argv[i], "--model_path")) model_path = argv[++i];
		else if (!strcmp(argv[i], "--hm_likelihood_online"))
        {
            ++i;
            if (!strcmp(argv[i], "yes") ||  !strcmp(argv[i], "YES"))
                harm_mean_likelihood_online = true;
        }
        else
        {
            printf("%s, unknown parameters, exit\n", argv[i]);
			char z;
			cout << ">>";
			cin >> z;
            exit(0);
        }
    }

	if (algorithm == NULL || data_path == NULL)
    {
        printf("Note that algorithm and data are not optional!\n");
        exit(0);
    }

    if (VERBOSE && (!strcmp(algorithm, "train") || (!strcmp(algorithm, "online")) ||
                    (!strcmp(algorithm, "train_dynamic")) || (!strcmp(algorithm, "online_dynamic"))))
    {
        printf("\nProgram starts with following parameters:\n");

        printf("algorithm:          = %s\n", algorithm);
        printf("data_path:          = %s\n", data_path);
		if (directory != NULL)
        printf("directory:          = %s\n", directory);

        printf("max_iter            = %d\n", max_iter);
        printf("random_seed         = %ld\n", seed);
        printf("vocab_size          = %d\n", m_vocab_size);
        printf("sampler_number      = %d\n", m_samplers);
        printf("gamma_a             = %.2f\n", gamma_a);
        printf("gamma_b             = %.2f\n", gamma_b);
        printf("alpha_a             = %.2f\n", alpha_a);
        printf("alpha_b             = %.2f\n", alpha_b);
        printf("eta                 = %.2f\n", eta);
        if (model_path != NULL)
        printf("saved model_path    = %s\n", model_path);
		if (harm_mean_likelihood_online)
			printf("computing likelihood with harmonic mean = yes\n");
		else
			printf("computing likelihood with harmonic mean = no\n");
    }

    // allocate the random number structure
    RANDOM_NUMBER = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(RANDOM_NUMBER, (long) seed); // init the seed


    if (!strcmp(algorithm, "train"))
    {
		run_train_batch_whole_process(m_samplers, m_vocab_size, gamma_a, gamma_b, alpha_a, 
								      alpha_b, eta, gamma, alpha, max_iter,
								      directory, data_path, model_path, 
									  save_topic_distribution);
    }

	if (!strcmp(algorithm, "online"))
	{
		run_online_process(m_samplers, m_vocab_size,
						   gamma_a, gamma_b, alpha_a,
						   alpha_b, eta, gamma, alpha, max_iter,
						   directory, data_path, data_path_online, 
						   model_path,
						   harm_mean_likelihood_online, save_topic_distribution,
						   m_online_documents);
	}
    
    if (!strcmp(algorithm, "train_dynamic"))
    {
        run_train_dynamic_batch_whole_process(m_samplers, m_vocab_size, gamma_a, gamma_b, alpha_a,
                                              alpha_b, eta, gamma, alpha, delta, max_iter,
                                              directory, data_path, model_path,
											  save_topic_distribution);
    }
    
    if (!strcmp(algorithm, "online_dynamic")) {
        run_online_dynamic_process(m_samplers, m_vocab_size,
                                   gamma_a, gamma_b, alpha_a,
                                   alpha_b, eta, gamma, alpha, delta, max_iter,
                                   directory, data_path, data_path_online,
                                   model_path,
								   harm_mean_likelihood_online,
								   save_topic_distribution, m_online_documents);
    }


    gsl_rng_free(RANDOM_NUMBER);
}
