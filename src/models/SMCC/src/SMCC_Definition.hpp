// SMCC definitions.

/**
 * @brief Secant bulk modulus definition:
 * 
 * \f[ K = \left(p^{\prime}/\Delta \epsilon_{e}^{vol}\right)  \left( \exp \left( \Delta \epsilon_{e}^{vol}/\kappa^{*} \right) -1 \right)    \f]
 * 
 * where \f$ p^{\prime} \f$ is the effective mean stress, \f$ \Delta \epsilon_{e}^{vol} \f$ is the elastic volumetric strain increment
 * and \f$ \kappa^{*} \f$ is the slope of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define SMCC_SECANT_BULK_MODULUS (p_prime/Delta_epsilon_e_vol)*(std::exp(Delta_epsilon_e_vol/kappa_star)-1)

/**
 * @brief Tangent bulk modulus definition:
 * 
 * \f[ K = \frac{p^{\prime}}{\kappa^{*}} \f]
 * 
 * where \f$ p^{\prime} \f$ is the mean effective stress and \f$ \kappa^{*} \f$ is the slope 
 * of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define SMCC_TANGENT_BULK_MODULUS p_prime/kappa_star

/**
 * @brief Shear modulus definition:
 * 
 * \f[ G = \frac{3 \left( 1-2\nu \right)K}{2 \left( 1+\nu \right)}\f]
 * 
 * where \f$ \nu \f$ is Poisson's ratio and \f$ K \f$ is the bulk modulus.
 */
#define SMCC_SHEAR_MODULUS (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu))

// Plastic definitions.

/**
 * @brief Yield surface function definition:
 * 
 * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c*s_ep \right)\f]
 * 
 * where \f$ q \f$ is the deviatoric stress, \f$ p^{\prime} \f$ is the mean effective stress, 
 * \f$ M\f$ is the frictional constant and \f$ p_c \f$ is the preconsolidation pressure.
 */
#define SMCC_YIELD std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c*s_ep)

/**
 * @brief Derivative of the yield surface with respect to the deviatoric stress.
 * 
 * \f[ \frac{\partial f}{\partial q} = 2q \f]
 * 
 * where \f$ q \f$ is the deviatoric stress.
 */
#define SMCC_DF_DQ 2*q

/**
 * @brief Derivatives of the plastic potential function with respect to the Lode angle:
 * 
 * \f[ \frac{\partial f}{\partial \theta} = 0 \f]
 */
#define SMCC_DF_DTHETA 0

/**
 * @brief Derivative of the yield surface with respect to the mean effective stress.
 * 
 * \f[ \frac{\partial f}{\partial p} = M^2\left(2 p^{\prime}-p_c s_{ep} \right) \f]
 * 
 * where \f$ M \f$ is the frictional constant, \f$ p_{\prime} \f$ is the mean effective stress, \f$ p_c \f$
 * is the pre-consolidation pressure and \f$ s_{ep} \f$ is the senstivity state variable.
 */
#define SMCC_DF_DP_PRIME std::pow(M,2)*(2*p_prime-p_c)

/**
 * @brief Derivatives of the plastic potential function with respect to the deviatoric stress:
 * 
 * \f[ \frac{\partial g}{\partial q} = 2q \f]
 */
#define SMCC_DG_DQ 2*q

/**
 * @brief Derivatives of the plastic potential function with respect to the Lode angle:
 * 
 * \f[ \frac{\partial g}{\partial \theta} = 0 \f]
 */
#define SMCC_DG_DTHETA 0

/**
 * @brief Derivative of the plastic potential function with respect to the effective mean stress:
 * 
 * \f[ \frac{\partial g}{\partial p} = M^2\left(2 p^{\prime}-p_c s_{ep} \right) \f]
 * 
 * where \f$ M \f$ is the frictional constant, \f$ p_{\prime} \f$ is the mean effective stress, \f$ p_c \f$
 * is the pre-consolidation pressure and \f$ s_{ep} \f$ is the senstivity state variable.
 */
#define SMCC_DG_DP_PRIME std::pow(M,2)*(2*p_prime-p_c)

/**
 * @brief Hardening modulus.
 *  
 * \f[ H = \frac{M^2 p^{\prime} p_c}{\lambda^*-\kappa^*} \operatorname{tr}\left( \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \right) \f]
 * 
 * where \f$ M \f$ is the frictional constant, f$ p_{\prime} \f$ is the mean effective stress, \f$ p_c \f$
 * is the pre-consolidation pressure, \f$ \lambda \f$ is the slope of the NCL, \f$ \kappa \f$ is the slope 
 * of the RCL and \f$ \frac{\partial f}{\partial \boldsymbol{\sigma}^{\prime}} \f$ are the derivatives of
 * the yield surface with respect to the effective stress state.
 */
#define SMCC_HARDENING_MODULUS ((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*s_ep*df_dsigma_prime.trace() + (std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2) + (A*2.0/3.0*(df_dsigma_prime* (df_dsigma_prime.transpose())).trace())))

/**
 * @brief State variable elastic update for void ratio:
 * 
 * \f[ e = e - e \Delta \epsilon_{vol, e}\f]
 */
#define SMCC_STATE_0_ELASTIC_UPDATE e-(e*Delta_epsilon_vol_e)

/**
 * @brief State variable elastic update for preconsolidation pressure:
 * 
 * \f[ p_{c} = constant \f]
 */
#define SMCC_STATE_1_ELASTIC_UPDATE p_c

/**
 * @brief State variable elastic update for sensitivity:
 * 
 * \f[ s_{ep} = constant \f]
 */
#define SMCC_STATE_2_ELASTIC_UPDATE s_ep

/**
 * @brief State variable plastic increment for void ratio:
 * 
 * \f[ \Delta e = -(1+e) \Delta \epsilon_{vol, p} \f]
 *
 */
#define SMCC_STATE_0_PLASTIC_INCREMENT -(1+e)*Delta_epsilon_vol_p

/**
 * @brief State variable plastic increment for preconsolidation pressure:
 * 
 * /f[ \Delta p_{c} = \Delta \lambda \frac{H_{p_{c}}}{M^2 p} /f]
 * 
 */
#define SMCC_STATE_1_PLASTIC_INCREMENT delta_lambda*(((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star))*s_ep*df_dsigma_prime.trace()/(std::pow(M,2)*p_prime*s_ep))

/**
 * @brief State variable plastic increment for sensitivity:
 * 
 * /f[ \Delta s_{ep} = \Delta \lambda \frac{H_{s_{ep}}}{M^2 p} /f]
 * 
 */
#define SMCC_STATE_2_PLASTIC_INCREMENT delta_lambda*(((std::pow(M,2)*p_prime*p_c)/(lambda_star-kappa_star)*-k*(s_ep-1.0)*std::sqrt((1-A)*std::pow(df_dsigma_prime.trace(),2) + (A*2.0/3.0*(df_dsigma_prime*(df_dsigma_prime.transpose())).trace())))/(std::pow(M,2)*p_prime*p_c))

