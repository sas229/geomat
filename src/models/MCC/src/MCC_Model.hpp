// Elastic definitions.

/**
 * @brief Secant bulk modulus definition:
 * 
 * \f[ K = \left(p^{\prime}/\Delta \epsilon_{e}^{vol}\right)  \left( \exp \left( \Delta \epsilon_{e}^{vol}/\kappa^{*} \right) -1 \right)    \f]
 * 
 * where \f$ p^{\prime} \f$ is the effective mean stress, \f$ \Delta \epsilon_{e}^{vol} \f$ is the elastic volumetric strain increment
 * and \f$ \kappa^{*} \f$ is the slope of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define BULK_MODULUS_SECANT K = (p_prime/delta_epsilon_e_vol)*(exp(delta_epsilon_e_vol/kappa_star)-1)

/**
 * @brief Tangent bulk modulus definition:
 * 
 * \f[ K = \frac{p^{\prime}}{\kappa^{*}} \f]
 * 
 * where \f$ p^{\prime} \f$ is the mean effective stress and \f$ \kappa^{*} \f$ is the slope 
 * of the recompression line in \f$ \ln \left( e \right)-\ln \left( p^{\prime} \right)\f$ space.
 */
#define BULK_MODULUS_TANGENT K = p_prime/kappa_star

/**
 * @brief Shear modulus definition:
 * 
 * \f[ G = \frac{3 \left( 1-2\nu \right)K}{2 \left( 1+\nu \right)}\f]
 * 
 * where \f$ \nu \f$ is Poisson's ratio and \f$ K \f$ is the bulk modulus.
 */
#define SHEAR_MODULUS G = (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu))

// Plastic definitions.

/**
 * @brief Yield surface function definition:
 * 
 * \f[ f = q^2 + M^2 p^{\prime}\left( p^{\prime}-p_c \right)\f]
 * 
 * where \f$ q \f$ is the deviatoric stress, \f$ p^{\prime} \f$ is the mean effective stress, 
 * \f$ M\f$ is the frictional constant and \f$ p_c \f$ is the preconsolidation pressure.
 */
#define YIELD f = pow(q,2) + pow(M,2)*p_prime*(p_prime-p_c)
