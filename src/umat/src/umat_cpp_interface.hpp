#ifndef UMAT_CPP_INTERFACE_H
#define UMAT_CPP_INTERFACE_H

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <memory>
#include <cassert>
#include <Eigen/Eigen>
#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include "Types.hpp"
#include "Model.hpp"
#include "models.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

/** @brief Interface to allow Abaqus to call umat material models via a C++ implementation.
 * 
 * This function allows Abaqus to call the library of user-defined materials defined within this library. 
 * Fortran passes variables by reference by default, hence pointers are passed for all variables. In the 
 * stress and strain arrays and in the matrices `ddsdde', `ddsddt', and `drplde', direct components are 
 * stored first, followed by shear components. There are `ndi' direct and `nshr' engineering shear components.
 * At a minimum the function must define `ddsdde', `stress', `statev'. 
 * 
 * The interface initializes Eigen objects using the pointers to arrays and passes them to the instantiated model
 * and on solution will remap the output to the input pointer arraysfor consumtpion by Abaqus. In this way, this
 * function simply acts as an interface between Abaqus and the generalised constitutive model objects derived from the 
 * Model base class provided by this framework. This approach allows the constitutive models to be easily wrapped 
 * for use in Python using pybind11.
 * 
 * @param[in,out] stress Stress tensor at the beginning of the increment to be updated.
 * @param[in,out] statev An array of solution dependent state variables to be updated.
 * @param[in,out] ddsdde Jacobian matrix to be updated.
 * @param[in,out] sse Specific elastic strain energy.
 * @param[in,out] spd Specific plastic dissipation.
 * @param[in,out] scd Specific creep dissipation.
 * @param[in,out] rpl Volumetric heat generation per unit time.
 * @param[in,out] ddsddt Variation of the stress increments with respect to the temperature.
 * @param[in,out] drplde Variation of rpl with respect to the strain increments.
 * @param[in,out] drpldt Variation of rpl with respect to the temperature.
 * @param[in,out] stran The total strains at the beginning of the increment.
 * @param[in,out] dstran Array of strain increments.
 * @param[in,out] time Value of step time (1) and total time (2) at the beginning of the increment.
 * @param[in,out] dtime Time increment.
 * @param[in,out] temp Temperature at the start of the increment.
 * @param[in,out] dtemp Temperature increment.
 * @param[in,out] predef Array of predefined field variables at this point at the start of the increment, interpolated from nodes.
 * @param[in,out] dpred Array of increments of predefined field variables.
 * @param[in,out] cmname User-defined material name (avoid ABQ_ prefix).
 * @param[in,out] ndi Number of direct stress components at this point.
 * @param[in,out] nshr Number of engineering shear stress components at this point.
 * @param[in,out] ntens Size of the stress or strain component array (`ndi' + `nshr').
 * @param[in,out] nstatv Number of solution-dependent state variables that are associated with this material type.
 * @param[in,out] props User-specified array of material constants associated with this user material.
 * @param[in,out] nprops Number of user-defined material constants associated with this user material.
 * @param[in,out] coords An array containing the coordinates of this point.
 * @param[in,out] drot Rotation increment matrix.
 * @param[in,out] pnewdt Ratio of suggested new time increment to the time increment being used.
 * @param[in,out] celent Characteristic element length, which is a typical length of a line across an element for a first-order element; it is half of the same typical length for a second-order element.
 * @param[in,out] dfgrd0 Array containing the deformation gradient at the beginning of the increment.
 * @param[in,out] dfgrd1 Array containing the deformation gradient at the end of the increment.
 * @param[in,out] noel Element number.
 * @param[in,out] npt Integration point number.
 * @param[in,out] layer Layer number.
 * @param[in,out] kspt Section point number within the current layer.
 * @param[in,out] kstep Procedure type key; nlgeom boolean; linear perturbation boolean.
 * @param[in,out] kinc Increment number.
 */
extern "C" void umat_cpp_interface(
    double *stress,
    double *statev,
    double *ddsdde,
    double *sse,
    double *spd,
    double *scd,
    double *rpl,
    double *ddsddt,
    double *drplde,
    double *drpldt,
    double *stran,
    double *dstran,
    double *time,
    double *dtime,
    double *temp,
    double *dtemp,
    double *predef,
    double *dpred,
    char *cmname,
    int *ndi,
    int *nshr,
    int *ntens,
    int *nstatv,
    double *props,
    int *nprops,
    double *coords,
    double *drot,
    double *pnewdt,
    double *celent,
    double *dfgrd0,
    double *dfgrd1,
    int *noel,
    int *npt,
    int *layer,
    int *kspt,
    int *kstep,
    int *kinc
    );

#ifdef __cplusplus
}
#endif

#endif