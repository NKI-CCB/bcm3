#include "Eigen/Dense"

bool initialize_likelihood(size_t num_variables, const char* const* variable_names)
{
    return true;
}

bool evaluate_log_probability(ptrdiff_t num_variables, const double* values, const char** variable_names, double* log_p)
{
    *log_p = 0.0;
    
    return true;
}
