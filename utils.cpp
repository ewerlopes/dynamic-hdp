#include "utils.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
using namespace std;

extern gsl_rng * RANDOM_NUMBER;

const double half_ln_2pi = 0.91893853320467267;

/**
 * given log(a) and log(b), return log(a + b)
 *
 */

double log_sum(double log_a, double log_b)
{
    double v;

    if (log_a < log_b)
        v = log_b+log(1 + exp(log_a-log_b));
    else
        v = log_a+log(1 + exp(log_b-log_a));
    return v;
}


// give a_1, ..., a_n,
// return log(exp(a_1)+...+exp(a_n))
double log_normalize(double * array, int nlen)
{
   const double log_max = 100.0; // the log(maximum in double precision), make sure it is large enough.
   int argmax;
   double max_val = max(array, nlen, &argmax); //get the maximum value in the array to avoid overflow
   double log_shift = log_max - log(nlen + 1.0) - max_val;
   double sum = 0.0;
   for (int i = 0; i < nlen; i++)
       sum += exp(array[i] + log_shift); //shift it

   double log_norm = log(sum) - log_shift;
   for (int i = 0; i < nlen; i++)
       array[i] -= log_norm; //shift it back

   return log_norm;
}

// the vector version
double log_normalize(vector<double> & vec, int nlen)
{
   const double log_max = 100.0; // the log(maximum in double precision), make sure it is large enough.
   int argmax;
   double max_val = max_vec(vec, nlen, &argmax); //get the maximum value in the array to avoid overflow
   double log_shift = log_max - log(nlen + 1.0) - max_val;
   double sum = 0.0;
   for (int i = 0; i < nlen; i++)
       sum += exp(vec[i] + log_shift); //shift it

   double log_norm = log(sum) - log_shift;
   for (int i = 0; i < nlen; i++)
       vec[i] -= log_norm; //shift it back

   return log_norm;
}

void log_normalize_non_zero_values(vector<double> & vec, int nlen)
{
    const double log_max = 100.0; // the log(maximum in double precision), make sure it is large enough.
    int argmax;
    vector<double> non_zero_values;
    non_zero_values.assign(nlen, 0);
    int non_zero_elements = 0;
    for (int i = 0; i < nlen; ++i) {
        if (abs(vec[i]) > 0) {
            non_zero_values[non_zero_elements] = vec[i];
            ++non_zero_elements;
        }
    }
    double max_val = max_vec(non_zero_values, non_zero_elements, &argmax); //get the maximum value in the array to avoid overflow
    double log_shift = log_max - log(nlen + 1.0) - max_val;
    double sum = 0.0;
    for (int i = 0; i < non_zero_elements; i++)
        sum += exp(non_zero_values[i] + log_shift); //shift it
    
    double log_norm = log(sum) - log_shift;
    for (int i = 0; i < nlen; i++) {
        if (abs(vec[i]) > 0) {
            vec[i] -= log_norm; //shift it back
        }
    }

    non_zero_values.clear();
}


double log_sum_exp_trick(const vector<double>& vector_of_logs)
// given log(x_i) return log(sum_i x_i) = log(sum_i exp(log(x_i)))
{
	double max_value = *max_element(vector_of_logs.begin(), vector_of_logs.end());

	double sum = 0;
	int vector_size = vector_of_logs.size();
	for (int i = 0; i < vector_size; ++i) {
		sum += (exp(vector_of_logs[i] - max_value));
	}
	double result = max_value + log(sum);

	return result;
}

/**
 * given log(a) and log(b), return log(a - b) a>b
 *
 */

double log_subtract(double log_a, double log_b)
{
    if (log_a < log_b) return -1000.0;

    double v;
    v = log_a + log(1 - exp(log_b-log_a));
    return v;
}

void invert_element_sign(vector<double>& input_vector)
{
	int vector_size = input_vector.size();
	for (int i = 0; i < vector_size; ++i) {
		input_vector[i] = -input_vector[i];
	}
}


/**
*
* check if file exisits
*/
//bool file_exists(const char * filename)
//{
//    if ( 0 == access(filename, R_OK))
//       return true;
//    return false;
//}

/**
*
* check if directory exisits
*/
bool dir_exists(const char * directory)
{
    struct stat st;
    if(stat(directory, &st) == 0)
        return true;
    return false;
}

/**
 * return factorial log((n-1+a)...(a))
 *
 **/

double log_factorial(int n, double a)
{
    if (n == 0) return 0.0;
    double v = lgamma(n+a) - lgamma(a);
    return v;
}

/**
 * return the cosine similarity
 *
 **/

double similarity(const int* v1, const int* v2, int n)
{
    double sim = 0.0, norm1 = 0.0, norm2 = 0.0;
    for (int i = 0; i < n; i ++)
    {
        sim += v1[i] * v2[i];
        norm1 += v1[i] * v1[i];
        norm2 += v2[i] * v2[i];
    }
    return sim/sqrt(norm1*norm2);
}


/// gsl_wrappers
double lgamma(double x) 
{
    return gsl_sf_lngamma(x);
}

unsigned int rmultinomial(const double* p, int n, double tot_p)
{
    int i;
    if (tot_p < 0)
    {
        tot_p = 0.0;
        for (i = 0; i < n; i ++) tot_p += p[i];
    }

    double u = runiform() * tot_p;
    double cum_p = 0.0;
    for (i = 0; i < n; i ++)
    {
        cum_p += p[i];
        if (u < cum_p) break;
    }
    return i;
}

double rgamma(double a, double b)
{
    return gsl_ran_gamma_mt(RANDOM_NUMBER, a, b);
}

double rbeta(double a, double b)
{
    return gsl_ran_beta(RANDOM_NUMBER, a, b);
}

unsigned int rbernoulli(double p)
{
    return gsl_ran_bernoulli(RANDOM_NUMBER, p);
}

double runiform()
{
    return gsl_rng_uniform_pos(RANDOM_NUMBER);
}

void rshuffle(void* base, size_t n, size_t size)
{
    gsl_ran_shuffle(RANDOM_NUMBER, base, n, size);
}

unsigned long int runiform_int(unsigned long int n)
{
    return gsl_rng_uniform_int(RANDOM_NUMBER, n);
}

void open_file(const char * filename, FILE * fileptr)
{
	printf("\nopening file %s\n", filename);

    fileptr = fopen(filename, "r");
}

void print_vector(const vector<double>& input_vector, const string& file_name)
{
    int vector_size = input_vector.size();
    if (vector_size == 0) {
        return;
    }
    
    ofstream out(file_name);
    out << input_vector[0];
    for (int i = 1; i < vector_size; ++i) {
        out << endl << input_vector[i];
    }
    out.close();
}

void print_matrix(const vector<int *>& input_matrix, int n_columns, const string& file_name)
{
    if ((input_matrix.size() == 0) || (n_columns == 0)) {
        return;
    }
    ofstream out(file_name);
    out << input_matrix[0][0];
    for (int d = 1; d < n_columns; ++d) {
        out << ", " << input_matrix[0][d];
    }
    for (int k = 1; k < input_matrix.size(); ++k) {
        out << endl;
        out << input_matrix[k][0];
        for (int d = 1; d < n_columns; ++d) {
            out << ", " << input_matrix[k][d];
        }
    }
    out.close();
}

void append_matrix(const vector<int>& input_vector, const string& file_name)
{
	if (input_vector.size() == 0) {
        return;
    }
    ofstream out;
	out.open(file_name, ofstream::out | ofstream::app);
    out << endl;
    out << input_vector[0];
    for (int k = 1; k < input_vector.size(); ++k) {
		out << ", " << input_vector[k];
    }
    out.close();
}

string add_postfix(const string& input_string, int postfix_num)
{
    string result = input_string;
    if (postfix_num != -1) {
        char* postfix = new char[5];
        postfix[0] = '_';
        postfix[1] = '0' + postfix_num / 100;
        postfix[2] = '0' + (postfix_num / 10) % 10;
        postfix[3] = '0' + postfix_num % 10;
        postfix[4] = '\0';

        result += postfix;
    }

    return result;
}

// end of the file
