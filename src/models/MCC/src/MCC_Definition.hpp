// Elastic definitions.

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
 * @brief Derivatives of the yield surface with respect to the effective stress state:
 * 
 * \f[ \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} = \frac{M^2\left(2 p^{\prime}-p_c \right)}{3} \mathbf{1} + 3 \operatorname{dev}(\boldsymbol{\sigma}^{\prime}) \f]
 */
#define DF_DSIGMA_PRIME pow(M,2)*(2*p_prime-p_c)*1.0/3.0*eye + 3.0*s

/**
 * @brief Derivatives of the plastic potential function with respect to the effective stress state. 
 * Equal to the derivates for the yield surface if associated flow:
 * 
 * \f[ \frac{\partial g}{\partial \boldsymbol{\sigma}^{\prime}} = \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \f]
 */
#define DG_DSIGMA_PRIME DF_DSIGMA_PRIME

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
 * @brief Void ratio elastic update:
 * 
 * \f[ e = e - e \delta \epsilon_{vol, e}\f]
 */
#define VOID_RATIO_ELASTIC_UPDATE e-(e*delta_epsilon_vol_e)

/**
 * @brief Preconsolidation pressure elastic update:
 * 
 * \f[ p_{c} = \text{constant}\f]
 */
#define PRECONSOLIDATION_ELASTIC_UPDATE p_c