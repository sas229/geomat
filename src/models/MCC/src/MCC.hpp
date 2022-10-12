#ifndef MCC_H
#define MCC_H

#include <string>
#include <iostream>
#include "Model.hpp"

class MCC : public Model {

    public: 

        /** @brief MCC model constructor. */
        MCC();
        /** @brief MCC model destructor. */
        virtual ~MCC() {}
        
        /** @brief Method to set the model parameters. */
        void set_parameters(std::vector<double> p);
        /** @brief Method to set state variables. */
        void set_state_variables(std::vector<double> s);

    protected:
    
        /** @brief Parameter: Poisson's ratio. */
        double nu;
        /** @brief Parameter: frictional constant. */
        double M;
        /** @brief State variable: intercept of NCL. */
        double N;
        /** @brief Parameter: slope of the NCL in ln(e)-ln(p') space. */
        double lambda_star;
        /** @brief Parameter: slope of the RCL in ln(e)-ln(p') space. */
        double kappa_star;

};

#endif 