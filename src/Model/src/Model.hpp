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

        /** @brief Virtual destructor for Model class. */
        virtual ~Model() {}

        // Setters.
        
        /** 
         * @brief Method to set Jacobian matrix. 
         */
        void set_jacobian(Eigen::MatrixXd j);

        /** 
         * @brief Virtual method to set model parameters. 
         */
        virtual void set_parameters(std::vector<double> s) {};

        /** 
         * @brief Method to set effective stress tensor in Voigt notation. 
         */
        void set_sigma_prime(Eigen::VectorXd s);

        /** 
         * @brief Virtual method to set state variables. 
         */
        virtual void set_state_variables(std::vector<double> s) {};

        /** 
         * @brief Method to set strain increment in Voigt notation. 
         */
        void set_strain_increment(Eigen::VectorXd s);
        
        // Getters.

        /**
         * @brief Method to get Jacobian matrix. 
         */
        Eigen::MatrixXd get_jacobian(void);

        /** 
         * @brief Method to get effective stress tensor in Voigt notation. 
         */
        Eigen::VectorXd get_sigma_prime(void);
        
        /** 
         * @brief Method to get state variables. 
         */
        std::vector<double> get_state_variables(void);

        /** 
         * @brief Method to get strain increment in Voigt notation.
         */
        Eigen::VectorXd get_strain_increment(void);

    // protected:

        // Setters.

        /** @brief Method to set name of model. */
        void set_name(std::string name);

        /** @brief Method to set number of model parameters. */
        void set_n_parameters(int i);
        
        /** @brief Method to set number of state variables. */
        void set_n_state_variables(int i);

        // Getters.

        /** @brief Method to get name of model. */
        std::string get_name(void);
        
        /** @brief Method to get number of model parameters. */
        int get_n_parameters(void);
        
        /** @brief Method to get number of state variables. */
        int get_n_state_variables(void);

        // Computers.

        /**
         * @brief Method to compute the cartesian stress tensor from principal stresses and directions as: 
         * \f[ \sigma_{i j}^{\prime} = T_{i j} S_{i j} T_{i j}^{T} \f] 
         * where \f$ T_{i j} \f$ is the principal stress direction tensor and \f$ S_{i j} \f$ is the principal stress tensor.
         * 
         * @param[in] T Principal stress direction tensor.
         * @param[in] S Principal stress tensor.
         * @returns Cartesian stress tensor.
         */
        Eigen::Matrix3d compute_cartesian_stresses(Eigen::Matrix3d T, Eigen::Matrix3d S);
        
        /** 
         * @brief Function to compute Lode's angle with cosine definition \f$ \theta_c \f$ as:
         * \f[ \theta_c=\arccos{\left[\frac{J_3}{2}\left(\frac{3}{J_2}\right)^{3/2}\right]} \f]
         * where \f$ J_2 \f$ and \f$ J_3 \f$ are the second and third deviatoric stress invariants, respectively.
         * This is converted to the Lode's angle for a sine definition \f$ \theta_s \f$ via:
         * \f[ \theta_s=\frac{\pi}{6}-\theta_c \f]
         * Similarly, the negative sine definition of Lode's angle \f$ \bar{\theta_s} \f$ is computed as:
         * \f[ \bar{\theta_s}=-\theta_s=-\frac{\pi}{6}+\theta_c\f]
         * 
         * @param[in] J_2 Second deviatoric stress invariant.
         * @param[in] J_3 Third deviatoric stress invariant.
         * @param[out] theta_c Lode's angle using cosine definition.
         * @param[out] theta_s Lode's angle using sine definition.
         * @param[out] theta_s_bar Lode's angle using negative sine definition.
         */
        void compute_lode(double J_2, double J_3, double &theta_c, double &theta_s, double &theta_s_bar);

        /** 
         * @brief Maximum shear stress calculated via:
         * \f[ \sigma_{max} = \frac{\max(\left|\sigma_1 - \sigma_2\right|, \left|\sigma_1 - \sigma_3\right|, \left|\sigma_2 - \sigma_3\right|)}{2}\f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the principal stresses. 
         * 
         * @param[in] sigma_1 Major principal stress.
         * @param[in] sigma_2 Intermediate principal stress.
         * @param[in] sigma_3 Minor principal stress.
         * @returns Maximum shear stress.
         */
        double compute_max_shear(double sigma_1, double sigma_2, double sigma_3);

        /** 
         * @brief von Mises stress calculated via:
         * \f[ \sigma_m = \sqrt{3J_{2}} \f]
         * where \f$ J_2 \f$ is the second deviatoric stress invariant.
         * 
         * @param[in] J_2 Second deviatoric stress invariant.
         * @returns von Mises stress.
         */
        double compute_mises_stress(double J_2);

        /** 
         * @brief Mean stress calculated via:
         * \f[ p=\frac{1}{3} \operatorname{tr}\left(\sigma_{i j}\right) \f]
         * where \f$\sigma_{i j} \f$ is the stress tensor.
         * 
         * @param[in] sigma Stress tensor.
         * @returns Mean stress.
         */
        double compute_p(Eigen::Matrix3d sigma);

        /** 
         * @brief Mean effective stress calculated via:
         * \f[ p^{\prime}=\frac{1}{3} \operatorname{tr}\left(\sigma_{i j}^{\prime}\right) \f]
         * where \f$\sigma_{i j}^{\prime} \f$ is the effective stress tensor.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @returns Mean effective stress.
         */
        double compute_p_prime(Eigen::Matrix3d sigma_prime);

        /** 
         * @brief Method to compute the principal stresses and directions. The principal stresses
         * are computed as the eigenvalues \f$ \lambda \f$ of the stress tensor via the characteristic
         * equation:
         * 
         * \f[ \left|\sigma_{i j}-\lambda \delta_{i j}\right|=-\lambda^3+I_1 \lambda^2-I_2 \lambda+I_3=0 \f]
         * 
         * The major principal stress \f$ \sigma_1 \f$ is the maximum eigenvalue, the intermediate principal stress  
         * \f$ \sigma_2 \f$ is the intermediate eigenvalue and the minor principal stress \f$ \sigma_3 \f$ is the 
         * is the minimum eigenvalue. These are subsequently assembled into the principal stress tensor #S and the 
         * associated eigenvectors are assembled into the principal stress direction tenor #T, which can be used to 
         * recover the cartesian stress state via compute_cartesian_stresses().
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @param[out] sigma_1 Major principal stress.
         * @param[out] sigma_2 Intermediate principal stress.
         * @param[out] sigma_3 Minor principal stress.
         * @param[out] S Principal stress tensor.
         * @param[out] T Principal stress direction tensor. 
         */
        void compute_principal_stresses(Eigen::Matrix3d sigma_prime, double &sigma_1, double &sigma_2, double &sigma_3, Eigen::Matrix3d &S, Eigen::Matrix3d &T);

        /** 
         * @brief Function to compute the deviatoric stress via:
         * \f[ q = \sqrt{\frac{1}{2} \left[
         * \left(\sigma_{1 1}^{\prime} -\sigma_{2 2}^{\prime}\right)^2 +
         * \left(\sigma_{2 2}^{\prime} -\sigma_{3 3}^{\prime}\right)^2 +
         * \left(\sigma_{3 3}^{\prime} -\sigma_{1 1}^{\prime}\right)^2 +
         * 6 \left( \tau_{1 2}^2 + \tau_{1 3}^2 + \tau_{2 3}^2 \right)
         * \right]} \f] 
         * where \f$ \sigma_{1 1}^\prime \f$, \f$ \sigma_{2 2}^\prime \f$, \f$ \sigma_{3 3}^\prime \f$, \f$ \tau_{1 2} \f$, \f$ \tau_{1 3} \f$
         * and \f$ \tau_{2 3} \f$ are the components of the effective stress tensor \f$ \sigma_{i j}^{\prime} \f$. 
         * 
         * @param[in] sigma
         * @returns Deviatoric stress. 
         */
        double compute_q(Eigen::Matrix3d sigma);

        /** 
         * @brief Function to compute the deviatoric stress tensor as: 
         * \f[ s = \sigma_{i j} - p \delta_{i j} \f]
         * where \f$ p \f$ is the mean stress and \f$ \delta_{i j} \f$ is the Kronecker delta.
         *  
         * @param[in] sigma Stress tensor.
         * @param[in] p Mean stress.
         * @returns Deviatoric stress tensor.
         */
        Eigen::Matrix3d compute_s(Eigen::Matrix3d sigma, double p);

        /** @brief Method to compute the stress invariants. The first, second and third stress invariants 
         * \f$ I_1 \f$, \f$ I_2 \f$ and \f$ I_3 \f$ are calculated as:
         * \f[ I_1=\operatorname{tr}\left(\sigma_{i j}\right) \f] 
         * \f[ I_2=\frac{1}{2}\left\{\left[\operatorname{tr}\left(\sigma_{i j}\right)\right]^2-\operatorname{tr}\left(\sigma_{i j}^2\right)\right\} \f]
         * \f[ I_3=\operatorname{det} \left(\sigma_{i j}\right) \f]
         * where \f$\sigma_{i j} \f$ is the total stress tensor. The first, second and third deviatoric stress invariants
         * \f$ J_1 \f$, \f$ J_2 \f$ and \f$ J_3 \f$ are calculated as:
         * \f[ J_1=\operatorname{tr}\left(s_{i j}\right)=0\f]
         * \f[ J_2=\frac{1}{2} \operatorname{tr}\left(s_{i j}^2\right) \f]
         * \f[ J_3=\operatorname{det} \left(s_{i j}\right) \f]
         * where \f$ s_{i j} \f$ is the deviatoric stress tensor, which can be computed via compute_s().
         * 
         * @param[in] sigma Total stress tensor.
         * @param[out] I_1 First stress invariant.
         * @param[out] I_2 Second stress invariant.
         * @param[out] I_3 Third stress invariant.
         * @param[out] J_1 First deviatoric stress invariant.
         * @param[out] J_2 Second deviatoric stress invariant.
         * @param[out] J_3 Third deviatoric stress invariant.
         */
        void compute_stress_invariants(Eigen::Matrix3d sigma, double &I_1, double &I_2, double &I_3, double &J_1, double &J_2, double &J_3);

        //  Updaters.
        
        /**
         * @brief Method to update the cartesian stress tensor class attribute #sigma_prime via compute_cartesian_stresses(). 
         */
        void update_cartesian_stresses(void);
        
        /** 
         * @brief Method to update Lode's angle class attributes #theta_c, #theta_s and #theta_s_bar via 
         * compute_lode(). 
         */
        void update_lode(void);

        /** 
         * @brief Method to update the maximum shear stress class attribute #max_shear via compute_max_shear(). 
         */
        void update_max_shear(void);

        /** 
         * @brief Method to update von Mises stress class attribute #mises_stress via compute_mises_stress(). 
         */
        void update_mises_stress(void);
        
        /** 
         * @brief Method to update the mean stress class attribute #p via compute_p(). 
         */
        void update_p(void);
        
        /** 
         * @brief Method to update the mean effective stress class attribute #p_prime via compute_p_prime(). 
         */
        void update_p_prime(void);

        /** 
         * @brief Method to compute the principal stresses #sigma_1, #sigma_2, #sigma_3, 
         * and principal stress and direction tensors #S and #T, 
         * via compute_principal_stresses(). 
         */
        void update_principal_stresses(void);

        /** 
         * @brief Method to update the deviatoric stress class attribute #q via compute_q(). 
         */
        void update_q(void);

        /** 
         * @brief Method to update the deviatoric stress tensor class attribute #s via compute_s(). 
         */
        void update_s(void);

        /** 
         * @brief Method to compute the stress invariant class attributes #I_1, #I_2, #I_3, #J_1, #J_2 and #J_3 via compute_stress_invariants(). 
         */
        void update_stress_invariants(void);

        // Members.

        /** 
         * @brief Name of model. 
         */
        std::string name;
        
        /** 
         * @brief Number of model parameters. 
         */
        int n_parameters;
        
        /** 
         * @brief Number of state variables. 
         */
        int n_state_variables;

        /** 
         * @brief Effective stress in Voigt notation:
         * \f[ \sigma^{\prime} = \left[\sigma_{1 1}^{\prime}, \sigma_{2 2}^{\prime}, \sigma_{3 3}^{\prime}, \tau_{1 2}, \tau_{1 3}, \tau_{2 3}\right]^T\f]. 
         */
        Eigen::VectorXd sigma_prime_v;

        /** 
         * @brief Effective stress tensor:
         * \f[ \sigma_{i j}^{\prime} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1 1}^{\prime} & \tau_{1 2} & \tau_{1 3} \\
         *      \tau_{2 1} & \sigma_{2 2}^{\prime} & \tau_{2 3} \\
         *      \tau_{3 1} & \tau_{3 2} & \sigma_{3 3}^{\prime}
         *      \end{array}\right] \f]
         * where \f$ \sigma_{1 1}^{\prime} \f$, \f$ \sigma_{2 2}^{\prime} \f$ and \f$ \sigma_{3 3}^{\prime} \f$ are the effective axial stresses, 
         * and  \f$ \tau_{1 2} \f$, \f$ \tau_{1 3} \f$ and \f$ \tau_{2 3} \f$ are the shear streses, respectively. 
         */
        Eigen::Matrix3d sigma_prime;

        /** 
         * @brief Total stress tensor:
         *  \f[ \sigma_{i j}=\sigma_{i j}^{\prime}+u\delta_{i j} \f]
         * where \f$ u \f$ is the pore pressure and \f$ \delta_{i j} \f$ is the Kronecker delta.
         */
        Eigen::Matrix3d sigma;

        /** 
         * @brief Deviatoric stress tensor. 
         */
        Eigen::Matrix3d s;

        /** 
         * @brief Strain increment in Voigt notation:
         * \f[ \Delta \epsilon = \left[\Delta \epsilon_{1 1}, \Delta \epsilon_{2 2}, \Delta \epsilon_{3 3}, \Delta \epsilon_{1 2}, \Delta \epsilon_{1 3}, \Delta \epsilon_{2 3}\right]^T\f]. 
         */
        Eigen::VectorXd delta_epsilon_v;

        /** 
         * @brief Strain increment tensor:
         * \f[ \Delta \epsilon_{i j} = 
         *      \left[\begin{array}{lll}
         *      \Delta \epsilon_{1 1} & \Delta \epsilon_{1 2} & \Delta \epsilon_{1 3} \\
         *      \Delta \epsilon_{2 1} & \Delta \epsilon_{2 2} & \Delta \epsilon_{2 3} \\
         *      \Delta \epsilon_{3 1} & \Delta \epsilon_{3 2} & \Delta \epsilon_{3 3}
         *      \end{array}\right] \f]
         * where \f$ \Delta \epsilon_{1 1} \f$, \f$ \Delta \epsilon_{2 2} \f$ and \f$ \Delta \epsilon_{3 3} \f$ are the axial strain increments, 
         * and  \f$ \Delta \epsilon_{1 2} \f$, \f$ \Delta \epsilon_{1 3} \f$ and \f$ \Delta \epsilon_{2 3} \f$ are the shear strain increments, respectively. 
         */
        Eigen::Matrix3d delta_epsilon;
        
        /** 
         * @brief Kronecker delta tensor:
         * \f[ \delta_{i j} = 
         *      \left[\begin{array}{lll}
         *      1 & 0 & 0 \\
         *      0 & 1 & 0 \\
         *      0 & 0 & 1
         *      \end{array}\right] \f] 
         */
        Eigen::Matrix3d delta {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
         
        /** 
         * @brief Pore pressure. Computed for undrained (bulk modulus) approach, otherwise provided as a state variable.
         */
        double u;

        /** 
         * @brief Jacobian matrix. 
         */
        Eigen::MatrixXd jacobian;

        /** 
         * @brief Array of state variables. 
         */
        std::vector<double> state;
        
        /** 
         * @brief First stress invariant.
         */
        double I_1; 

        /** 
         * @brief Second stress invariant. 
         */
        double I_2; 

        /** 
         * @brief Third stress invariant.
         */
        double I_3; 

        /** 
         * @brief First deviatoric stress invariant.
         */
        double J_1; 

        /** 
         * @brief Second deviatoric stress invariant.
         */
        double J_2; 

        /** 
         * @brief Third deviatoric stress invariant.
         */
        double J_3; 

        /** 
         * @brief Mean stress calculated by compute_p().
         */
        double p; 

        /** 
         * @brief Effective mean stress calculated by compute_p_prime().
         */
        double p_prime; 

        /**
         * @brief Deviatoric stress calculated by compute_q(). 
         */
        double q; 

        /** 
         * @brief Principal stress tensor \f$ S_{i j} \f$: 
         * \f[ S_{i j} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1} & 0 & 0 \\
         *      0 & \sigma_{2} & 0 \\
         *      0 & 0 & \sigma_{3}
         *      \end{array}\right] \f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the major, intermediate and minor principal stresses, respectively. 
         */
        Eigen::Matrix3d S;

        /** 
         * @brief Major principal stress. 
         */
        double sigma_1; 

        /** 
         * @brief Intermediate principal stress. 
         */
        double sigma_2; 

        /** 
         * @brief Minor principal stress. 
         */
        double sigma_3;

        /** 
         * @brief Principal stress directions tensor \f$ T_{i j} \f$. 
         */
        Eigen::Matrix3d T;

        /** 
         * @brief Mises stress:
         * \f[ \sigma_{m} = \sqrt{3 J_2}\f]
         * where \f$ J_2 \f$ is the second deviatoric stress invariant. 
         */
        double mises_stress;

        /** 
         * @brief Maximum shear stress:
         * \f[ \sigma_{max} = \frac{\max(\left|\sigma_1 - \sigma_2\right|, \left|\sigma_1 - \sigma_3\right|, \left|\sigma_2 - \sigma_3\right|)}{2}\f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the principal stresses. 
         */
        double max_shear;  

        /** 
         * @brief Constant \f$ \pi \f$. 
         */
        double pi = 2*std::acos(0.0);

        /** 
         * @brief Lode angle from cosine definition. 
         */
        double theta_c;

        /** 
         * @brief Lode angle from sine definition. 
         */
        double theta_s;

        /** 
         * @brief Lode angle from negative sine definition. 
         */
        double theta_s_bar;
};

#endif