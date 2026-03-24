#ifndef UDSM_CPP_INTERFACE_H
#define UDSM_CPP_INTERFACE_H

#include <memory>
#include <mutex>
#include <cstring>
#include <Eigen/Eigen>
#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include "Types.hpp"
#include "Model.hpp"
#include "ModelLookup.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

/** @brief Interface to allow Plaxis to call UDSM material models via a C++ implementation.
 *
 * This function implements the Plaxis User-Defined Soil Model (UDSM) interface, allowing
 * constitutive models defined within this framework to be used in Plaxis analyses. Fortran
 * passes variables by reference by default, hence pointers are passed for all variables.
 *
 * The interface initialises Eigen objects using the pointers to arrays and passes them to
 * the instantiated model. On solution, the output is remapped to the input pointer arrays
 * for consumption by Plaxis. The sign convention is converted from the Plaxis
 * compression-negative convention to the compression-positive convention used internally.
 *
 * @param[in] IDTask Identification of the task: 1 = initialise state variables,
 *            2 = calculate constitutive stresses, 3 = create effective stress-related
 *            constitutive matrix, 4 = return number of state variables, 5 = inquire
 *            matrix properties, 6 = create elastic material matrix.
 * @param[in] iMod Model number (1-based index into the list of available models).
 * @param[in] IsUndr Flag indicating drained (0) or undrained (1) analysis.
 * @param[in] iStep Current calculation step number.
 * @param[in] iTer Current iteration number within the step.
 * @param[in] iEl Element number.
 * @param[in] Int Integration point number.
 * @param[in] X X-coordinate of the integration point.
 * @param[in] Y Y-coordinate of the integration point.
 * @param[in] Z Z-coordinate of the integration point.
 * @param[in] Time0 Time at the start of the current step.
 * @param[in] dTime Time increment for the current step.
 * @param[in] Props Array of material parameter values.
 * @param[in] Sig0 Stress tensor at the start of the increment (6 components).
 * @param[in] Swp0 Excess pore water pressure at the start of the increment.
 * @param[in] StVar0 State variables at the start of the increment.
 * @param[in] dEps Strain increment tensor (6 components).
 * @param[in,out] D Constitutive matrix (6x6, stored column-major).
 * @param[in,out] Bulk_W Bulk modulus of water (for undrained behaviour).
 * @param[in,out] Sig Updated stress tensor (6 components).
 * @param[in,out] Swp Updated excess pore water pressure.
 * @param[in,out] StVar Updated state variables.
 * @param[in,out] ipl Plasticity indicator (0 = elastic, 1 = plastic).
 * @param[in,out] nStat Number of state variables used by the model.
 * @param[in,out] NonSym Flag indicating non-symmetric constitutive matrix (0 = symmetric, 1 = non-symmetric).
 * @param[in,out] iStrsDep Flag indicating stress-dependent stiffness (0 = no, 1 = yes).
 * @param[in,out] iTimeDep Flag indicating time-dependent behaviour (0 = no, 1 = yes).
 * @param[in,out] iTang Flag controlling tangent stiffness calculation.
 * @param[in] iPrjDir Project directory path string.
 * @param[in] iPrjLen Length of the project directory path string.
 * @param[in,out] iAbort Abort flag (set to 1 to signal an error to Plaxis).
 */
extern "C" void udsm_cpp_interface(
    int *IDTask,
    int *iMod,
    int *IsUndr,
    int *iStep,
    int *iTer,
    int *iEl,
    int *Int,
    double *X,
    double *Y,
    double *Z,
    double *Time0,
    double *dTime,
    double *Props,
    double *Sig0,
    double *Swp0,
    double *StVar0,
    double *dEps,
    double *D,
    double *Bulk_W,
    double *Sig,
    double *Swp,
    double *StVar,
    int *ipl,
    int *nStat,
    int *NonSym,
    int *iStrsDep,
    int *iTimeDep,
    int *iTang,
    char *iPrjDir,
    int *iPrjLen,
    int *iAbort
);

#ifdef __cplusplus
}
#endif

#endif
