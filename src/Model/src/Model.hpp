#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector> 
#include <string>
#include <plog/Log.h>
#include <Eigen/Eigen>

/** @brief The Model class contains methods and attributes that are common to all genera of constitutive model. 
 * The Model class is the base class from which all constitutive models are derived. */
class Model {

    public:

        /** @brief Method to set stress tensor in Voigt notation. */
        void set_sigma(Eigen::VectorXd s);
        
        /** @brief Method to set Jacobian matrix. */
        void set_jacobian(Eigen::MatrixXd j);
        
        /** @brief Method to get stress tensor in Voigt notation. */
        Eigen::VectorXd get_sigma(void);
        
        /** @brief Method to get Jacobian matrix. */
        Eigen::MatrixXd get_jacobian(void);
        
        /** @brief Method to get state variables. */
        std::vector<double> get_state_variables(void);

        /** @brief Virtual method to set state variables. */
        virtual void set_state_variables(std::vector<double> s) {};
        
        /** @brief Virtual method to set model parameters. */
        virtual void set_parameters(std::vector<double> s) {};

        /** @brief Virtual destructor for Model class. */
        virtual ~Model() {}
        
    protected:

        /** @brief Method to set name of model. */
        void set_name(std::string name);

        /** @brief Method to set number of model parameters. */
        void set_n_parameters(int i);
        
        /** @brief Method to set number of state variables. */
        void set_n_state_variables(int i);

        /** @brief Method to get name of model. */
        std::string get_name(void);
        
        /** @brief Method to get number of model parameters. */
        int get_n_parameters(void);
        
        /** @brief Method to get number of state variables. */
        int get_n_state_variables(void);

        /** @brief Method to compute the stress invariants. */
        void compute_stress_invariants(void);

        /** @brief Method to compute the principal stresses and directions. */
        void compute_principal_stresses(void);
        
        /** @brief Method to compute the cartesian stress tensor from principal stresses and directions, derived via: 
         * \f[ \sigma_{i j}^{\prime} = T_{i j} S_{i j} T_{i j}^{T} \f] 
         * where \f$ T_{i j} \f$ is the principal stress direction tensor and \f$ S_{i j} \f$ is the principal stress tensor.*/
        void compute_cartesian_stresses(void);

        /** @brief Name of model. */
        std::string name;
        
        /** @brief Number of model parameters. */
        int n_parameters;
        
        /** @brief Number of state variables. */
        int n_state_variables;

        /** @brief Effective stress in Voigt notation:
         * \f[\left[\sigma_{11}^{\prime},\sigma_{22}^{\prime},\sigma_{33}^{\prime},\tau_{12},\tau_{13},\tau_{23}\right]^T\f]. */
        Eigen::VectorXd sigma_prime_v;

        /** @brief Effective stress tensor:
         * \f[ \sigma_{i j}^{\prime} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1 1}^{\prime} & \tau_{1 2} & \tau_{1 3} \\
         *      \tau_{2 1} & \sigma_{2 2}^{\prime} & \tau_{2 3} \\
         *      \tau_{3 1} & \tau_{3 2} & \sigma_{3 3}^{\prime}
         *      \end{array}\right] \f]
         * where \f$ \sigma_{1 1}^{\prime} \f$, \f$ \sigma_{2 2}^{\prime} \f$ and \f$ \sigma_{3 3}^{\prime} \f$ are the effective axial stresses, 
         * and  \f$ \tau_{1 2} \f$, \f$ \tau_{1 3} \f$ and \f$ \tau_{2 3} \f$ are the shear streses, respectively. */
        Eigen::Matrix<double, 3, 3> sigma_prime;

        /** @brief Total stress tensor:
         *  \f[ \sigma_{i j}=\sigma_{i j}^{\prime}+u\delta_{i j} \f]
         * where \f$ u \f$ is the pore pressure and \f$ \delta_{i j} \f$ is the Kronecker delta.*/
        Eigen::Matrix<double, 3, 3> sigma;

        /** @brief Deviatoric stress tensor. 
         * \f[ s = \sigma_{i j} - p \f]
         * where \f$ p \f$ is the mean stress.
        */
        Eigen::Matrix<double, 3, 3> s;
        
        /** @brief Kronecker delta tensor:
         * \f[ \delta_{i j} = 
         *      \left[\begin{array}{lll}
         *      1 & 0 & 0 \\
         *      0 & 1 & 0 \\
         *      0 & 0 & 1
         *      \end{array}\right] \f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the major, intermediate and minor principal stresses, respectively. */
        Eigen::Matrix<double, 3, 3> delta {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
         
        /** @brief Pore pressure. Computed for undrained (bulk modulus) approach, otherwise provided as a state variable.*/
        double u;

        /** @brief Jacobian matrix. */
        Eigen::MatrixXd jacobian;
        /** @brief Array of state variables. */
        std::vector<double> state;
        
        /** @brief First stress invariant calculated via: 
         * \f[ I_1=\operatorname{tr}\left(\sigma_{i j}\right) \f] 
         * where \f$\sigma_{i j} \f$ is the stress tensor.*/
        double I_1; 

        /** @brief Second stress invariant calculated via: 
         * \f[ I_2=\frac{1}{2}\left\{\left[\operatorname{tr}\left(\sigma_{i j}\right)\right]^2-\operatorname{tr}\left(\sigma_{i j}^2\right)\right\} \f]. */
        double I_2; 

        /** @brief Third stress invariant calculated via:
         * \f[ I_3=\operatorname{det} \left(\sigma_{i j}\right) \f]
         * where \f$\sigma_{i j} \f$ is the stress tensor.*/
        double I_3; 

        /** @brief First deviatoric stress invariant calculated via:
         * \f[ J_1=\operatorname{tr}\left(s_{i j}\right)=0\f]
         * where \f$s_{i j}\f$ is the deviatoric stress tensor.
         */
        double J_1; 

        /** @brief Second deviatoric stress invariant calculated via:
         * \f[ J_2=\frac{1}{2} \operatorname{tr}\left(s_{i j}^2\right) \f]
         * where \f$s_{i j}\f$ is the deviatoric stress tensor.
         */
        double J_2; 

        /** @brief Third deviatoric stress invariant calculated via:
         * \f[ J_3=\operatorname{det} \left(s_{i j}\right) \f]
         * where \f$s_{i j} \f$ is the deviatoric stress tensor.*/
        double J_3; 

        /** @brief Mean stress calculated via:
         * \f[ p=\frac{1}{3} \operatorname{tr}\left(\sigma_{i j}\right) \f]
         * where \f$\sigma_{i j} \f$ is the stress tensor.*/
        double p; 

        /** @brief Effective mean stress calculated via:
         * \f[ p^{\prime}=\frac{1}{3} \operatorname{tr}\left(\sigma_{i j}^{\prime}\right) = p - u \f]
         * where \f$\sigma_{i j}^{\prime} \f$ is the effective stress tensor and \f$ u \f$ is the pore pressure.*/
        double p_prime; 

        /** @brief Deviatoric stress calculated via: 
         * \f[ q=\sqrt{\frac{1}{2}\left[\left(\sigma_1-\sigma_2\right)^2+\left(\sigma_2-\sigma_3\right)^2+\left(\sigma_1-\sigma_3\right)^2\right]} \f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the major, intermediate and minor principal stresses, respectively.
        */
        double q; 

        /** @brief Principal stress tensor: 
         * \f[ S_{i j} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1} & 0 & 0 \\
         *      0 & \sigma_{2} & 0 \\
         *      0 & 0 & \sigma_{3}
         *      \end{array}\right] \f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the major, intermediate and minor principal stresses, respectively. */
        Eigen::Matrix3d S;

        /** @brief Major principal stress. */
        double sigma_1; 

        /** @brief Intermediate principal stress. */
        double sigma_2; 

        /** @brief Minor principal stress. */
        double sigma_3;

        /** @brief Principal stress directions tensor \f$ T_{i j} \f$. */
        Eigen::Matrix3d T;

        /** @brief Mises stress:
         * \f[ \sigma_{m} = \sqrt{3 J_2}\f]
         * where \f$ J_2 \f$ is the second deviatoric stress invariant. */
        double mises;

        /** @brief Maximum shear stress:
         * \f[ \sigma_{max} = \frac{\max(\left|\sigma_1 - \sigma_2\right|, \left|\sigma_1 - \sigma_3\right|, \left|\sigma_2 - \sigma_3\right|)}{2}\f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the principal stresses. */
        double max_shear;  

        /** @brief Constant \f$ \pi \f$. */
        double pi = 2*std::acos(0.0);

        /** @brief Lode angle from cosine relationship:
         * \f[ \theta_c=\arccos{\left[\frac{J_3}{2}\left(\frac{3}{J_2}\right)^{3/2}\right]} \f]
         * where \f$ J_2 \f$ and \f$ J_3 \f$ are the second and third deviatoric stress invariants, respectively. */
        double theta_c;

        /** @brief Lode angle from sine relationship:
         * \f[ \theta_s=\frac{\pi}{6}-\theta_c \f]
         * where \f$ \theta_c \f$ is the Lode angle from the cosine relationship. */
        double theta_s;

        /** @brief Lode angle from negative sine relationship:
         * \f[ \bar{\theta_s}=-\theta_s=-\frac{\pi}{6}+\theta_c\f]
         * where \f$ \theta_s \f$ and \f$ \theta_c \f$ are the Lode angle derived from the sine and cosine relationships, respectively. */
        double theta_s_bar;
};

#endif