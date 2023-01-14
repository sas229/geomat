#ifndef CHECKS_H
#define CHECKS_H

#include <plog/Log.h>
#include <string>
#include <math.h>

namespace Checks {
    /**
     * @brief Function to check model inputs.
     * 
     * @param name Model name.
     * @param parameters_size Size of parameters vector.
     * @param state_size Size of state variables vector.
     * @param parameters_required Number of parameters required.
     * @param state_required Number of state variables required.
     */
    void check_inputs(std::string name, int parameters_size, int state_size, int parameters_required, int state_required);

    /**
     * @brief Check the elastic parameters.
     * 
     * @param K Bulk modulus.
     * @param G Secant modulus.
     */
    void check_elastic_parameters(double K, double G);
}

#endif