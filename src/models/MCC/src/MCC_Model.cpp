// Elastic definition.
#define SECANT_BULK_MODULUS (p_prime/delta_epsilon_vol)*(std::exp(delta_epsilon_vol/kappa_star)-1)
#define TANGENT_BULK_MODULUS p_prime/kappa_star
#define SHEAR_MODULUS (3.0*(1-2.0*nu)*K)/(2.0*(1.0+nu))

// Yield surface definition.
#define YIELD std::pow(q,2) + std::pow(M,2)*p_prime*(p_prime-p_c)