#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <plog/Log.h>
#include <Eigen/Eigen>

class Model {
    public:
        /** @brief Method to set stress tensor. */
        void set_stress(Eigen::VectorXd stress);

        /** @brief Method to get stress tensor. */
        Eigen::VectorXd get_stress(void);
        
    protected:

        /** @brief Method to set name of model. */
        void set_name(std::string name);
        /** @brief Method to set number of model parameters. */
        void set_nparams(int nparams);
        /** @brief Method to set number of state variables. */
        void set_nstatev(int nstatev);

        /** @brief Method to get name of model. */
        std::string get_name(void);
        /** @brief Method to get number of model parameters. */
        int get_nparams(void);
        /** @brief Method to get number of state variables. */
        int get_nstatev(void);

        /** @brief Name of model. */
        std::string _name;
        /** @brief Number of model parameters. */
        int _nparams;
        /** @brief Number of state variables. */
        int _nstatev;
        /** @brief Number of variables in stress tensor. */
        int _ntens;
        /** @brief Stress tensor in Voigt notation \f$(\sigma_{11},\sigma_{22},\sigma_{33},\sigma_{23},\sigma_{13},\sigma_{12})\f$. */
        Eigen::VectorXd _stress;
};

#endif