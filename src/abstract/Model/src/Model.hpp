#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector> 
#include <string>
#include <plog/Log.h>
#include <Eigen/Eigen>
#include "Types.hpp"

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
        void set_jacobian(Jacobian j);

        /** 
         * @brief Method to set effective stress tensor in Voigt notation. 
         */
        void set_sigma_prime(Voigt sigma_prime_tilde);

        /** 
         * @brief Method to set strain increment in Voigt notation. 
         */
        void set_strain_increment(Voigt delta_epsilon);
        
        // Getters.

        /**
         * @brief Method to get Jacobian matrix. 
         */
        Jacobian get_jacobian(void);

        /** 
         * @brief Method to get effective stress tensor in Voigt notation. 
         */
        Voigt get_sigma_prime(void);

        /**
         * @brief Method to solve stress increment given the current strain increment.
         * 
         * @note Must be overriden by constitutive behaviour classes (e.g. Elastic, Elastoplastic).
         */
        virtual void solve(void) {};

        double get_p_prime(void);

        /** 
         * @brief Mean effective stress calculated via:
         * \f[ p^{\prime}=\frac{1}{3} \operatorname{tr}\left( \boldsymbol{\sigma}^{\prime} \right) \f]
         * where \f$ \boldsymbol{\sigma}^{\prime} \f$ is the effective stress tensor.
         * 
         * @param[in] sigma_prime Effective stress tensor.
         * @returns Mean effective stress.
         */
        double compute_p_prime(Cauchy sigma_prime);

        Cauchy compute_dq_dsigma_prime(Cauchy sigma_prime);

    protected:

        // Setters.

        /** @brief Method to set name of model. */
        void set_name(std::string name);

        // Getters.

        /** @brief Method to get name of model. */
        std::string get_name(void);

        // Computers.

        /**
         * @brief Method to compute the cartesian stress tensor from principal stresses and directions as: 
         * \f[ \boldsymbol{\sigma}^{\prime} = \mathbf{T S T}^\intercal \f] 
         * where \f$ \mathbf{T} \f$ is the principal stress direction tensor and \f$ \mathbf{S} \f$ is the principal stress tensor.
         * 
         * @param[in] T Principal stress direction tensor.
         * @param[in] S Principal stress tensor.
         * @returns Cartesian stress tensor.
         */
        Cauchy compute_cartesian_stresses(Cauchy T, Cauchy S);

        /** 
         * @brief Method to compute the volumetric strain increment:
         * \f[ \Delta \epsilon_{vol} = \operatorname{tr} \left( \boldsymbol{\Delta \epsilon} \right)\f] 
         * where \f$ \boldsymbol{\Delta \epsilon} \f$ is the strain increment tensor. 
         * 
         * @param[in] delta_epsilon Strain increment tensor.
         * @return Volumetric strain increment.
         */
        double compute_delta_epsilon_vol(Cauchy delta_epsilon);

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
         * \f[ p=\frac{1}{3} \operatorname{tr}\left( \boldsymbol{\sigma} \right) \f]
         * where \f$ \boldsymbol{\sigma} \f$ is the total stress tensor.
         * 
         * @param[in] sigma Total stress tensor.
         * @returns Mean stress.
         */
        double compute_p(Cauchy sigma);

        /** 
         * @brief Method to compute the principal stresses and directions. The principal stresses
         * are computed as the eigenvalues \f$ \lambda \f$ of the stress tensor via the characteristic
         * equation:
         * 
         * \f[ -\lambda^3+I_1 \lambda^2-I_2 \lambda+I_3=0 \f]
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
         * @param[out] R Principal stress direction tensor. 
         * @param[out] S Principal stress tensor.
         */
        void compute_principal_stresses(Cauchy sigma_prime, double &sigma_1, double &sigma_2, double &sigma_3, Cauchy &R, Cauchy &S);

        /** 
         * @brief Function to compute the deviatoric stress via:
         * \f[ q = \sqrt{\frac{1}{2} \left[
         * \left(\sigma_{1 1} -\sigma_{2 2}\right)^2 +
         * \left(\sigma_{2 2} -\sigma_{3 3}\right)^2 +
         * \left(\sigma_{3 3} -\sigma_{1 1}\right)^2 +
         * 6 \left( \tau_{1 2}^2 + \tau_{1 3}^2 + \tau_{2 3}^2 \right)
         * \right]} \f] 
         * where \f$ \sigma_{1 1} \f$, \f$ \sigma_{2 2} \f$, \f$ \sigma_{3 3} \f$, \f$ \tau_{1 2} \f$, \f$ \tau_{1 3} \f$
         * and \f$ \tau_{2 3} \f$ are the components of the stress tensor \f$ \boldsymbol{\sigma} \f$. 
         * 
         * @param[in] sigma
         * @returns Deviatoric stress. 
         */
        double compute_q(Cauchy sigma);

        /** 
         * @brief Function to compute the deviatoric stress tensor as: 
         * \f[ s =  \boldsymbol{\sigma} - p \mathbf{I} \f]
         * where \f$ p \f$ is the mean stress and \f$ \mathbf{I} \f$ is the identity matrix.
         *  
         * @param[in] sigma Stress tensor.
         * @param[in] p Mean stress.
         * @returns Deviatoric stress tensor.
         */
        Cauchy compute_s(Cauchy sigma, double p);

        /** 
         * @brief Function to compute the total stress tensor as: 
         * \f[ \boldsymbol{\sigma} = \boldsymbol{\sigma}^{\prime} + u \mathbf{I} \f]
         * where \f$ \boldsymbol{\sigma}^{\prime} \f$ is the effective stress tensor, 
         * \f$ u \f$ is the pore pressure and \f$ \mathbf{I} \f$ is the identity matrix.
         *  
         * @param[in] sigma_prime Stress tensor.
         * @param[in] u Pore pressure.
         * @returns Total stress tensor.
         */
        Cauchy compute_sigma(Cauchy sigma_prime, double u);

        /** @brief Method to compute the stress invariants. The first, second and third stress invariants 
         * \f$ I_1 \f$, \f$ I_2 \f$ and \f$ I_3 \f$ are calculated as:
         * \f[ I_1=\operatorname{tr}\left( \boldsymbol{\sigma}\right) \f] 
         * \f[ I_2=\frac{1}{2}\left\{\left[\operatorname{tr}\left( \boldsymbol{\sigma}\right)\right]^2-\operatorname{tr}\left( \boldsymbol{\sigma}^2\right)\right\} \f]
         * \f[ I_3=\operatorname{det} \left( \boldsymbol{\sigma}\right) \f]
         * where \f$ \boldsymbol{\sigma} \f$ is the total stress tensor. The first, second and third deviatoric stress invariants
         * \f$ J_1 \f$, \f$ J_2 \f$ and \f$ J_3 \f$ are calculated as:
         * \f[ J_1=\operatorname{tr}\left(\mathbf{s}\right)=0\f]
         * \f[ J_2=\frac{1}{2} \operatorname{tr}\left(\mathbf{s}^2\right) \f]
         * \f[ J_3=\operatorname{det} \left(\mathbf{s}\right) \f]
         * where \f$ \mathbf{s} \f$ is the deviatoric stress tensor, which can be computed via compute_s().
         * 
         * @param[in] sigma Total stress tensor.
         * @param[out] I_1 First stress invariant.
         * @param[out] I_2 Second stress invariant.
         * @param[out] I_3 Third stress invariant.
         * @param[out] J_1 First deviatoric stress invariant.
         * @param[out] J_2 Second deviatoric stress invariant.
         * @param[out] J_3 Third deviatoric stress invariant.
         */
        void compute_stress_invariants(Cauchy sigma, double &I_1, double &I_2, double &I_3, double &J_1, double &J_2, double &J_3);

        //  Updaters.
        
        /**
         * @brief Method to update the cartesian stress tensor class attribute #sigma_prime via compute_cartesian_stresses(). 
         */
        void update_cartesian_stresses(void);

        /**
         * @brief Method to update the volumetric strain increment tensor class attribute #delta_epsilon_vol via compute_delta_epsilon_vol(). 
         */
        void update_delta_epsilon_vol(void);
        
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
         * @brief Method to update the total stress tensor class attribute #s via compute_sigma(). 
         */
        void update_sigma(void);

        /** 
         * @brief Method to compute the stress invariant class attributes #I_1, #I_2, #I_3, #J_1, #J_2 and #J_3 via compute_stress_invariants(). 
         */
        void update_stress_invariants(void);

        Voigt to_voigt(Cauchy cauchy);

        Cauchy to_cauchy(Voigt voigt);

        // Members.

        /** 
         * @brief Name of model. 
         */
        std::string name;

        /** 
         * @brief Effective stress in Voigt notation:
         * \f[ \boldsymbol{\tilde{\sigma}}^{\prime} = \left[\sigma_{1 1}^{\prime}, \sigma_{2 2}^{\prime}, \sigma_{3 3}^{\prime}, \tau_{1 2}, \tau_{1 3}, \tau_{2 3}\right]^T\f]
         */
        Voigt sigma_prime_tilde;

        /** 
         * @brief Effective stress tensor:
         * \f[  \boldsymbol{\sigma}^{\prime} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1 1}^{\prime} & \tau_{1 2} & \tau_{1 3} \\
         *      \tau_{2 1} & \sigma_{2 2}^{\prime} & \tau_{2 3} \\
         *      \tau_{3 1} & \tau_{3 2} & \sigma_{3 3}^{\prime}
         *      \end{array}\right] \f]
         * where \f$ \sigma_{1 1}^{\prime} \f$, \f$ \sigma_{2 2}^{\prime} \f$ and \f$ \sigma_{3 3}^{\prime} \f$ are the effective axial stresses, 
         * and  \f$ \tau_{1 2} \f$, \f$ \tau_{1 3} \f$ and \f$ \tau_{2 3} \f$ are the shear streses, respectively. 
         */
        Cauchy sigma_prime = Cauchy::Zero();

        /** 
         * @brief Total stress tensor:
         *  \f[  \boldsymbol{\sigma}= \boldsymbol{\sigma}^{\prime}+u\mathbf{I} \f]
         * where \f$ u \f$ is the pore pressure and \f$ \mathbf{I} \f$ is the identity matrix.
         */
        Cauchy sigma = Cauchy::Zero();

        /** 
         * @brief Deviatoric stress tensor. 
         */
        Cauchy s = Cauchy::Zero();

        /** 
         * @brief Strain increment in Voigt notation:
         * \f[ \Delta \tilde{\epsilon} = \left[\Delta \epsilon_{1 1}, \Delta \epsilon_{2 2}, \Delta \epsilon_{3 3}, \Delta \epsilon_{1 2}, \Delta \epsilon_{1 3}, \Delta \epsilon_{2 3}\right]^T\f]. 
         */
        Voigt delta_epsilon_tilde;

        /** 
         * @brief Strain increment tensor:
         * \f[ \Delta \boldsymbol{\epsilon} = 
         *      \left[\begin{array}{lll}
         *      \Delta \epsilon_{1 1} & \Delta \epsilon_{1 2} & \Delta \epsilon_{1 3} \\
         *      \Delta \epsilon_{2 1} & \Delta \epsilon_{2 2} & \Delta \epsilon_{2 3} \\
         *      \Delta \epsilon_{3 1} & \Delta \epsilon_{3 2} & \Delta \epsilon_{3 3}
         *      \end{array}\right] \f]
         * where \f$ \Delta \epsilon_{1 1} \f$, \f$ \Delta \epsilon_{2 2} \f$ and \f$ \Delta \epsilon_{3 3} \f$ are the axial strain increments, 
         * and  \f$ \Delta \epsilon_{1 2} \f$, \f$ \Delta \epsilon_{1 3} \f$ and \f$ \Delta \epsilon_{2 3} \f$ are the shear strain increments, respectively. 
         */
        Cauchy delta_epsilon = Cauchy::Zero();

        /** 
         * @brief Volumetric strain increment:
         * \f[ \Delta \epsilon_{vol} = \operatorname{tr} \left(\Delta \boldsymbol{\epsilon} \right)\f] 
         * where \f$ \Delta \boldsymbol{\epsilon} \f$ is the strain increment tensor. 
         */
        double delta_epsilon_vol;
        
        /** 
         * @brief Identity matrix:
         * \f[ \mathbf{I} = 
         *      \left[\begin{array}{lll}
         *      1 & 0 & 0 \\
         *      0 & 1 & 0 \\
         *      0 & 0 & 1
         *      \end{array}\right] \f] 
         */
        Cauchy eye = Cauchy::Identity();
         
        /** 
         * @brief Pore pressure. Computed for undrained (bulk modulus) approach, otherwise provided as a state variable.
         */
        double u;

        /** 
         * @brief Jacobian matrix. 
         */
        Jacobian jacobian = Jacobian::Zero();
    
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
         * @brief Principal stress tensor \f$ \mathbf{S} \f$: 
         * \f[ \mathbf{S} = 
         *      \left[\begin{array}{lll}
         *      \sigma_{1} & 0 & 0 \\
         *      0 & \sigma_{2} & 0 \\
         *      0 & 0 & \sigma_{3}
         *      \end{array}\right] \f]
         * where \f$ \sigma_1 \f$, \f$ \sigma_2 \f$ and \f$ \sigma_3 \f$ are the major, intermediate and minor principal stresses, respectively. 
         */
        Cauchy S = Cauchy::Zero();

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
         * @brief Principal stress directions tensor \f$ R_{i j} \f$. 
         */
        Cauchy R = Cauchy::Zero();

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

        /**
         * @brief Boolean indicating whether the current strain increment has been solved.
         */
        bool solved;

        /**
         * @brief Derivative of the mean stress with respect to the effective stress state:
         * \f[ \frac{\partial p}{\partial \boldsymbol{\sigma}^{\prime}} = \frac{1}{3} \boldsymbol{I} \f]
         * where \f$ \boldsymbol{I} \f$ is the identity matrix. 
         */
        Cauchy dp_dsigma_prime = 1.0/3.0*eye;

        /**
         * @brief Derivative of the deviatoric stress with respect to the effective stress state:
         * \f[ \frac{\partial q}{\partial \boldsymbol{ \sigma^{\prime}}} = \frac{3}{2q} \left[\begin{array}{lll}
         *      s_{11} & 2\tau_{12} & 2\tau_{13} \\
         *      2\tau_{21} & s_{22} & 2\tau_{23} \\
         *      2\tau_{31} & 2\tau_{32} & s_{33}
         *      \end{array}\right] \f]
         * where \f$ s_{ii} \f$ are the deviatoric stress components, \f$ \tau_{ii} \f$ are the shear stresses and \f$ q \f$ is the deviatoric stress.
         */
        Cauchy dq_dsigma_prime;
};

#endif