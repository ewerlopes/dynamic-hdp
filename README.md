# dynamic Hierarchical Dirichlet Process
___

This is the C++ implementation of a dynamic Hierarchical Dirichlet process topic model. For the vanilla batch HDP, implementation by Chong Wang and David Blei is used.

The model is introduced in:

O.Isupova, D.Kuzin, L.Mihaylova. **Dynamic Hierarchical Dirichlet Process for Abnormal Behaviour Detection in Video**. Accepted to the 19th International Conference on Information Fusion 2016 and ICML 2016 Anomaly Detection Workshop.

Citation:
```
@inproceedings{isupova2016dynamic,
  title={Dynamic Hierarchical Dirichlet Process for abnormal behaviour detection in video},
  author={Isupova, Olga and Kuzin, Danil and Mihaylova, Lyudmila},
  booktitle={Information Fusion (FUSION), 2016 19th International Conference on},
  pages={750--757},
  year={2016},
  organization={IEEE}
}
```

Written by Olga Isupova and Danil Kuzin (`{o.isupova, dkuzin1}@sheffield.ac.uk`)

---

### DEPENDENCIES

* GNU Scientific Library.

    For ubuntu:

```
sudo apt-get install libgsl-dev
```
---

### BUILDING

After cloning the repository, on terminal issue the commands:

```
mkdir build && cd build
cmake .. && make
```
---

### DATA FORMAT:
a file where each line is of the form (the LDA-C format):

```
[M] [term_1]:[count] [term_2]:[count] ...  [term_N]:[count]
```
where [M] is the number of unique terms in the document, and the
[count] associated with each term is how many times that term appeared
in the document.

---

### RUNNING PARAMETERS:

1. Vanilla HDP for batch learning:

```
./HDP_experiments --algorithm train --data data_file --directory dir_to_output --model_path save_model_for_later_use --sampler_number number_of_posterior_samples --gamma first_level_concentration_parameter --alpha second_level_concentration_parameter --delta weight_for_global_topic_usage --eta Dirichlet_parameter
```

This runs the vanilla HDP on the given train data with the given number of posterior samples. The models are saved for the later use for online inference and anomaly detection.


2. Vanilla HDP for online learning:
```
./HDP_experiments --algorithm online --data data_file --online_data testing_data_file --directory dir_to_output --model_path save_model_for_later_use --sampler_number number_of_posterior_samples --hm_likelihood_online yes --gamma first_level_concentration_parameter --alpha second_level_concentration_parameter --delta weight_for_global_topic_usage --eta Dirichlet_parameter
```

This runs the vanilla HDP for online inference on the given testing data. It also computes the likelihood for each testing document and saves it in the "predictive_hm_likelihood" file  in the `dir_to_output` directory.

3. Dynamic HDP for batch learning:
```
./HDP_experiments --algorithm train_dynamic --data data_file --directory dir_to_output --model_path save_model_for_later_use --sampler_number number_of_posterior_samples --gamma first_level_concentration_parameter --alpha second_level_concentration_parameter --delta weight_for_global_topic_usage --eta Dirichlet_parameter
```

This runs the dynamic HDP on the given train data with the given number of posterior samples. The models are saved for the later use for online inference and anomaly detection.

4. Dynamic HDP for online learning:
```./HDP_experiments --algorithm online_dynamic --data data_file --online_data testing_data_file --directory dir_to_output --model_path save_model_for_later_use --sampler_number number_of_posterior_samples --hm_likelihood_online yes --gamma first_level_concentration_parameter --alpha second_level_concentration_parameter --delta weight_for_global_topic_usage --eta Dirichlet_parameter
```

This runs the dynamic HDP for online inference on the given testing data. It also computes the likelihood for each testing document and saves it in the "predictive_hm_likelihood" file  in the `dir_to_output` directory.

---

### SYNTHETIC DATA

In the paper the synthetic data is used for testing the algorithms. The train synthetic data is saved in the file `data/synthetic_dynamic_train_LDA_format.dat`, the testing synthetic data is saved in the file `data/synthetic_dynamic_test_LDA_format.dat`.

In order to get the results from the paper it is required:

1. Train batch vanilla HDP:
```
./HDP_experiments --algorithm train --data synthetic_dynamic_train_LDA_format.dat --directory synthetic_vanilla_HDP --model_path synthetic_vanilla_HDP_batch_model --vocab_size 25 --sampler_number 5 --random_seed 2000 --gamma 2 --delta 3 --eta 0.2 --alpha 1.5
```

2. Run online vanilla HDP:
```
./HDP_experiments --algorithm online --data synthetic_dynamic_train_LDA_format.dat --directory synthetic_vanilla_hdp --online_data synthetic_dynamic_test_LDA_format.dat --model_path synthetic_vanilla_hdp_batch_model --vocab_size 25 --sampler_number 5 --hm_likelihood_online yes --random_seed 2000 --online_data_size 1000 --alpha 1.5 --gamma 2 --delta 0.5 --eta 0.2
```

3. Train batch dynamic HDP:
```
./HDP_experiments --algorithm train_dynamic --data synthetic_dynamic_train_LDA_format.dat --directory synthetic_dynamic_HDP --model_path synthetic_dynamic_batch_model --vocab_size 25 --sampler_number 5 --random_seed 2000 --gamma 2 --delta 0.5 --eta 0.2 --alpha 1.5
```
4. Run online dynamic HDP:
```
./HDP_experiments --algorithm online_dynamic --data synthetic_dynamic_train_LDA_format.dat --directory synthetic_dynamic_HDP --online_data synthetic_dynamic_test_LDA_format.dat --model_path synthetic_dynamic_batch_model --vocab_size 25 --sampler_number 5 --hm_likelihood_online yes --random_seed 2000 --online_data_size 1000 --alpha 1.5 --gamma 2 --delta 0.5 --eta 0.2
```

5. Run matlab script `script/handle_results.m` for graphical representation

### Copyright
(C) Copyright 2010, Chong Wang and David Blei. Written by [Chong Wang](http://www.cs.princeton.edu/~chongw/index.html).

(C) Copyright 2016, Olga Isupova and Danil Kuzin.
