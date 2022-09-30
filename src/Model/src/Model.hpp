#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <plog/Log.h>
#include <Eigen/Eigen>

/** @brief Vector with six indexes. */
typedef Eigen::Vector<double, 6> Vector6d;

class Model {
    public:
    
        // Setters.
        
        /** @brief Method to set name of model. */
        void set_name(std::string name);
        /** @brief Method to set number of model parameters. */
        void set_nparams(int nparams);
        /** @brief Method to set stress tensor. */
        void set_stress(Vector6d stress);
        
        // Getters.

        /** @brief Method to get name of model. */
        std::string get_name(void);
        /** @brief Method to get number of model parameters. */
        int get_nparams(void);
        /** @brief Method to get stress tensor. */
        Vector6d get_stress(void);
        
    private:
        /** @brief Name of model. */
        std::string _name;
        /** @brief Number of model parameters. */
        int _nparams;
        /** @brief Stress tensor in Voigt notation \f$(\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{12},\sigma_{23},\sigma_{31})\f$. */
        Vector6d _stress;
};

#endif