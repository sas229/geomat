
#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Core>

/** 
 * @brief Constant \f$ \pi \f$. 
 */
const double pi = 2*std::acos(0.0);

/** 
 * @brief Eigen vector with six indices for storing the Voigt tensor. 
 */
typedef Eigen::Vector<double, 6> Voigt;

/** 
 * @brief Eigen square matrix with three indices for storing the Cauchy tensor. 
 */
typedef Eigen::Matrix<double, 3, 3> Cauchy;

/** 
 * @brief Eigen square matrix with six indices for storing the constitutive matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Constitutive;

/** 
 * @brief Eigen square matrix with six indices for storing the Jacobian matrix. 
 */
typedef Eigen::Matrix<double, 6, 6> Jacobian;

/**
 * @brief Vector of parameters.
 */
typedef Eigen::VectorXd Parameters;

/**
 * @brief Vector of state variables.
 */
typedef Eigen::VectorXd State;

/**
 * @brief Vector of hardening moduli.
 */
typedef Eigen::VectorXd HardeningModuli;

/**
 * @brief Vector of state variable factors.
 */
typedef Eigen::VectorXd StateFactors;

/**
 * @brief Structure containing the model derivatives.
 */
struct Derivatives {
    Cauchy df_dsigma_prime; /** @brief Derivatives of the yield function with respect to the stress state.*/
    Cauchy dg_dsigma_prime; /** @brief Derivatives of the plastic potential function with respect to the stress state.*/
    HardeningModuli H_s;    /** @brief Hardening moduli vector.*/
    StateFactors B_s;       /** @brief State factor vector.*/

    // Setters.
    void set_df_dsigma_prime(Cauchy df_dsigma_prime) {this->df_dsigma_prime = df_dsigma_prime;};
    void set_dg_dsigma_prime(Cauchy dg_dsigma_prime) {this->dg_dsigma_prime = dg_dsigma_prime;};
    void set_H_s(HardeningModuli H_s) {this->H_s = H_s;};
    void set_B_s(StateFactors B_s) {this->B_s = B_s;};

    // Getters.
    Cauchy get_df_dsigma_prime(void) {return this->df_dsigma_prime;};
    Cauchy get_dg_dsigma_prime(void) {return this->dg_dsigma_prime;};
    HardeningModuli get_H_s(void) {return this->H_s;};
    StateFactors get_B_s(void) {return this->B_s;};
};

/**
 * @brief Yield function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam state_f State variables.
 */
typedef std::function<double(Cauchy sigma_prime_f, State state_f)> YieldFunction;

/**
 * @brief Trial stress function binding.
 * 
 * @tparam sigma_prime Effective stress tensor.
 * @tparam state State variables.
 * @tparam delta_epsilon_tilde Strain increment in Voigt form.
 * @tparam sigma_prime_n Trial effective stress tensor.
 * @tparam state_n Trial state variables.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Voigt delta_epsilon_tilde, Cauchy &sigma_prime_n, State &state_n)> TrialFunction;

/**
 * @brief Constitutive matrix function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor.
 * @tparam Delta_epsilon Strain increment in Cuachy form.
 */
typedef std::function<Constitutive(Cauchy sigma_prime_f, Cauchy Delta_epsilon)> ConstitutiveMatrixFunction;

/**
 * @brief Derivative computation function binding.
 * 
 * @tparam sigma_prime Effective stress tensor.
 * @tparam state State variables.
 */
typedef std::function<Derivatives(Cauchy sigma_prime, State state)> DerivativeFunction;

/**
 * @brief Plastic increment function binding.
 * 
 * @tparam sigma_prime_f Effective stress tensor in Cauchy form.
 * @tparam state_f State variables.
 * @tparam Delta_epsilon_tilde_p_dT Plastic strain increment in Voigt form.
 * @tparam Delta_sigma_prime Increment in the effective stress state in Voigt form.
 * @tparam Delta_state Increment in the state variables.
 */
typedef std::function<void(Cauchy sigma_prime, State state, Voigt Delta_epsilon_p_dT, Voigt &Delta_sigma_prime, State &Delta_state)> PlasticIncrementFunction;

/**
 * @brief Yield surface drift function binding.
 * 
 * @tparam sigma_prime_u Uncorrected effective stress tensor in Cauchy form.
 * @tparam state_u Uncorrected state variables.
 * @tparam sigma_prime_c Corrected effective stress tensor in Cauchy form.
 * @tparam state_c Corrected state variables.
 * @tparam ITS_YSC Number of yield surface correction iterations.
 */
typedef std::function<void(Cauchy sigma_prime_u, State state_u, Cauchy &sigma_prime_c, State &state_c, int &ITS_YSC)> CheckDriftFunction;


/**
 * @brief Object containing functions to be bound to model specific implementation.
 */
struct ModelFunctions {
    YieldFunction compute_f;
    TrialFunction compute_trial_increment;
    ConstitutiveMatrixFunction compute_D_e;
    DerivativeFunction compute_derivatives;
    PlasticIncrementFunction compute_plastic_increment;
    CheckDriftFunction check_yield_surface_drift;
};

/**
 * @brief Stress integration settings with defaults.
 */
struct Settings {
    std::string solver = "Explicit";
    std::string method = "ModifiedEuler";
    double FTOL = 1e-8;     /** @brief Yield surface tolerance. */
    double LTOL = 1e-6;     /** @brief Unload-reload tolerance. */
    double STOL = 1e-4;     /** @brief Stress integration tolerance. */
    double EPS = 1e-16;     /** @brief Numerical precision tolerance. */
    double DT_MIN = 1e-6;   /** @brief Minimum pseudo-time increment for stress integration procedure. */
    int MAXITS_YSI = 10;    /** @brief Maximum number of Pegasus method iterations. */
    int MAXITS_YSC = 30;    /** @brief Maximum number of yield surface correction iterations. */
    int NSUB = 10;          /** @brief Maximum number of unload-reload intersection substeps.*/

    // Setters.
    void set_FTOL(double FTOL) {this->FTOL = FTOL;};
    void set_LTOL(double LTOL) {this->LTOL = LTOL;};
    void set_STOL(double STOL) {this->STOL = STOL;};
    void set_EPS(double EPS) {this->EPS = EPS;};
    void set_DT_MIN(double DT_MIN) {this->DT_MIN = DT_MIN;};
    void set_MAXITS_YSI(int MAXITS_YSI) {this->MAXITS_YSI = MAXITS_YSI;};
    void set_MAXITS_YSC(int MAXITS_YSC) {this->MAXITS_YSC = MAXITS_YSC;};
    void set_NSUB(int NSUB) {this->NSUB = NSUB;};

    // Getters.
    double get_FTOL(void) {return FTOL;};
    double get_LTOL(void) {return LTOL;};
    double get_STOL(void) {return STOL;};
    double get_EPS(void) {return EPS;};
    double get_DT_MIN(void) {return DT_MIN;};
    int get_MAXITS_YSI(void) {return MAXITS_YSI;};
    int get_MAXITS_YSC(void) {return MAXITS_YSC;};
    int get_NSUB(void) {return NSUB;};
};

#endif


         
         
         