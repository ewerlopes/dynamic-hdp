%% load true answers
load('true_answers_for_synthetic.mat');

%% load "true" roc curve
load('roc_curve_for_synthetic.mat');

%% load algorithm answers
dynamic_hm_likelihood = csvread('./synthetic_dynamic_HDP/predictive_hm_likelihood');

vanilla_hm_likelihood = csvread('./synthetic_vanilla_hdp/predictive_hm_likelihood');


%% auc
[x_dyn_hm, y_dyn_hm, t, auc_dyn_hm] = perfcurve(true_answers, dynamic_hm_likelihood, 0);
[x_van_hm, y_van_hm, t, auc_van_hm] = perfcurve(true_answers, vanilla_hm_likelihood, 0);

%% plotting
plot(x_dyn_hm, y_dyn_hm, 'b');
hold on
plot(x_van_hm, y_van_hm, '--r')
plot(x_true, y_true, ':m');
xlabel('False positive rate');
ylabel('True positive rate');
title('ROC-curves on synthetic data');
legend('Dynamic HDP', 'HDP', '"True" model')




