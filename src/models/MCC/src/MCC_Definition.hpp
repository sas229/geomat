// MCC definitions.

/**
 * @brief Secant bulk modulus definition:
 * 
 * \f[ K = \left(p^{\prime}/\Delta \epsilon_{e}^{vol}\right)  \left( \exp \left( \Delta \epsilon_{e}^{vol}/\kappa^{*} \right) -1 \right)    \f]
 * 
 * where \f$ p^{\prime} \f$ is the effective mean stress, \f$ \Delta \epsilon_{e}^{vol} \f$ is the elastic volumetric strain increment
 * and \f$ \kappa^{*} \f$ is the slope of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define BULK_MODULUS_SECANT (p_prime/delta_epsilon_e_vol)*(exp(delta_epsilon_e_vol/kappa_star)-1)

/**
 * @brief Tangent bulk modulus definition:
 * 
 * \f[ K = \frac{p^{\prime}}{\kappa^{*}} \f]
 * 
 * where \f$ p^{\prime} \f$ is the mean effective stress and \f$ \kappa^{*} \f$ is the slope 
 * of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define BULK_MODULUS_TANGENT p_prime/kappa_star

/**
 * @brief Shear modulus definition:
 * 
 * \f[ G = \frac{3 \left( 1-2\nu \right)K}{2 \left( 1+\nu \right)}\f]
 * 
 * where \f$ \nu \f$ is Poisson's ratio and \f$ K \f$ is the bulk modulus.
 */
#define SHEAR_MODULUS (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu))

// Plastic definitions.

/**
 * @brief Yield surface function definition:
 * 
 * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c \right)\f]
 * 
 * where \f$ q \f$ is the deviatoric stress, \f$ p^{\prime} \f$ is the mean effective stress, 
 * \f$ M\f$ is the frictional constant and \f$ p_c \f$ is the preconsolidation pressure.
 */
#define YIELD pow(q,2) + pow(M,2)*p_prime*(p_prime-p_c)

/**
 * @brief Derivatives of the yield surface with respect to the deviatoric stress:
 * 
 * \f[ \frac{\partial f}{\partial q} = 2q \f]
 */
#define DF_DQ 2*q

/**
 * @brief Derivatives of the yield surface with respect to the mean effective stress:
 * 
 * \f[ \frac{\partial f}{\partial p} = M^2\left(2 p^{\prime}-p_c \right) \f]
 */
#define DF_DP_PRIME pow(M,2)*(2*p_prime-p_c)

/**
 * @brief Derivatives of the plastic potential function with respect to the deviatoric stress:
 * 
 * \f[ \frac{\partial g}{\partial q} = 2q \f]
 */
#define DG_DQ 2*q

/**
 * @brief Derivative of the plastic potential function with respect to the effective mean stress:
 * 
 * \f[ \frac{\partial g}{\partial p^{\prime}} = M^2 \left( 2 p^{\prime} p_c \right) \f]
 */
#define DG_DP_PRIME pow(M,2)*(2*p_prime-p_c)

/**
 * @brief Hardening modulus.
 *  
 * \f[ H = \frac{M^2 p^{\prime} p_c}{\lambda^*-\kappa^*} \operatorname{tr}\left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \right) \f]
 */
#define HARDENING_MODULUS (pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*df_dsigma_prime.trace()

/**
 * @brief State variable elastic update for void ratio:
 * 
 * \f[ e = e - e \delta \epsilon_{vol, e}\f]
 */
#define STATE_0_ELASTIC_UPDATE e-(e*delta_epsilon_vol_e)

/**
 * @brief State variable elastic update for preconsolidation pressure:
 * 
 * \f[ p_{c} = constant \f]
 */
#define STATE_1_ELASTIC_UPDATE p_c

/**
 * @brief State variable plastic increment for void ratio:
 * 
 * \f[ \Delta e = -(1+e) \delta \epsilon_{vol, p} \f]
 *
 */
#define STATE_0_PLASTIC_INCREMENT -(1+e)*delta_epsilon_vol_p

/**
 * @brief State variable plastic increment for preconsolidation pressure:
 * 
 * /f[ \Delta p_{c} = \Delta \lambda \frac{H}{M^2 p} /f]
 * 
 */
#define STATE_1_PLASTIC_INCREMENT delta_lambda*H/(pow(M,2)*p_prime)