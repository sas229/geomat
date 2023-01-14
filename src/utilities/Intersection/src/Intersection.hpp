#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <plog/Log.h>
#include "Types.hpp"

class Intersection {
    
    public:

        virtual double compute_f_alpha(Cauchy sigma_prime, State state, double alpha, Voigt Delta_epsilon_tilde) = 0;

        virtual bool check_unload_reload(Cauchy sigma_prime, State state) = 0;

        double compute_alpha(Cauchy sigma_prime, State state, Voigt Delta_epsilon_tilde, double FTOL);

};

#endif