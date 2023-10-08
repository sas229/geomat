#ifndef LinearElastic_H
#define LinearElastic_H

#include "Elastic.hpp"
#include "Tensor.hpp"
#include "Types.hpp"

class LinearElastic : public Elastic {

    public: 

        /** 
         * @brief LinearElastic model constructor. 
         * 
         * @param[in] parameters Vector of parameters.
         * @param[in] state Vector of state variables.
         * @param[in] log_severity Severity of message to log.
         */
        LinearElastic(Parameters parameters, State state, std::string log_severity="none");

        /** 
         * @brief LinearElastic model destructor. 
         */
        virtual ~LinearElastic() {}

        /**
         * @brief Get the state variables vector.
         * 
         * @return State
         */
        State get_state_variables(void) {return state;};

        void set_state_variables(State new_state) {this->state = new_state;};

    protected:

        /**
         * @brief Overriden method to compute the elastic matrix.
         * 
         * @param sigma_prime Effective stress tensor.
         * @param Delta_epsilon Strain increment.
         * @return Constitutive 
         * 
         * Isotropic linear elastic constitutive matrix:
         * 
         * \f[ D_e = \left[\begin{array}{cccccc}
                     K + \frac{4}{3}G & K - \frac{2}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
                     K - \frac{2}{3}G & K + \frac{4}{3}G & K - \frac{2}{3}G & 0 & 0 & 0 \\
                     K - \frac{2}{3}G & K - \frac{2}{3}G & K + \frac{4}{3}G & 0 & 0 & 0 \\
                     0 & 0 & 0 & G & 0 & 0 \\
                     0 & 0 & 0 & 0 & G & 0 \\
                     0 & 0 & 0 & 0 & 0 & G
                     \end{array}\right] \f]
         * where \f$ K \f$ is the bulk modulus and \f$ G \f$ is the shear modulus.
         */
        Constitutive compute_D_e(Cauchy sigma_prime, Cauchy Delta_epsilon=Cauchy::Zero()) override;

        /** 
         * @brief Parameters. 
         */
        Parameters parameters;

        /** 
         * @brief State variables. 
         */
        State state;

        /** 
         * @brief Parameter: bulk modulus. 
         */
        const double &K = parameters[0];
        
        /** 
         * @brief Parameter: shear modulus. 
         */
        const double &G = parameters[1];

        /**
         * @brief Number of required parameters.
         */
        int parameters_required = 2;

        /**
         * @brief Number of required state variables.
         */
        int state_required = 0;

};

#endif 
