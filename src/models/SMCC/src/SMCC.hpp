#ifndef SMCC_H
#define SMCC_H

#include <string>
#include <iostream>
#include "Model.hpp"

class SMCC : public Model{

    public: 

        /** 
         * @brief SMCC model constructor. 
         */
        SMCC();

        /** 
         * @brief SMCC model destructor. 
         */
        virtual ~SMCC() {}
        
    protected:

        /** 
         * @brief State variable: lambda_star. 
         */
        double lambda_star;
        
};

#endif
