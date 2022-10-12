#ifndef SMCC_H
#define SMCC_H

#include <string>
#include <iostream>
#include "Model.hpp"

class SMCC : public Model{

    public: 

        /** @brief SMCC model constructor. */
        SMCC();
        /** @brief SMCC model destructor. */
        virtual ~SMCC() {}

        /** @brief Method to set state variables. */
        void set_state_variables(std::vector<double> s);
        
    protected:

        /** @brief State variable: lambda_star. */
        double lambda_star;
        
};

#endif
