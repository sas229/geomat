#ifndef UMAT_H
#define UMAT_H

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <memory>

#include <Eigen/Eigen>

#include <g3log/g3log.hpp>
#include <g3log/logworker.hpp>

#include "Model.hpp"
#include "MCC.hpp"
#include "SMCC.hpp"

/** Abaqus umat interface.
 * 
 * This function allows Abaqus to call the library of user-defined materials defined within this library. 
 * Fortran passes variables by reference by default, hence pointers are passed for all variables. In the 
 * stress and strain arrays and in the matrices `ddsdde', `ddsddt', and `drplde', direct components are 
 * stored first, followed by shear components. There are `ndi' direct and `nshr' engineering shear components.
 * At a minimum the function must define `ddsdde', `stress', `statev'. 
 * 
 * @brief Interface to allow Abaqus to call umat material models via a C++ implementation.
 * @param stress Stress tensor at the beginning of the increment to be updated.
 * @param statev An array of solution dependent state variables to be updated.
 * @param ddsdde Jacobian matrix to be updated.
 * @param sse Specific elastic strain energy.
 * @param spd Specific plastic dissipation.
 * @param scd Specific creep dissipation.
 * @param rpl Volumetric heat generation per unit time.
 * @param ddsddt Variation of the stress increments with respect to the temperature.
 * @param drplde Variation of rpl with respect to the strain increments.
 * @param drpldt Variation of rpl with respect to the temperature.
 * @param stran The total strains at the beginning of the increment.
 * @param dstran Array of strain increments.
 * @param time Value of step time (1) and total time (2) at the beginning of the increment.
 * @param dtime Time increment.
 * @param temp Temperature at the start of the increment.
 * @param predef Array of predefined field variables at this point at the start of the increment, interpolated from nodes.
 * @param dpred Array of increments of predefined field variables.
 * @param cmname User-defined material name (avoid ABQ_ prefix).
 * @param ndi Number of direct stress components at this point.
 * @param nshr Number of engineering shear stress components at this point.
 * @param ntens Size of the stress or strain component array (`ndi' + `nshr').
 * @param nstatv Number of solution-dependent state variables that are associated with this material type.
 * @param props User-specified array of material constants associated with this user material.
 * @param nprops Number of user-defined material constants associated with this user material.
 * @param coords An array containing the coordinates of this point.
 * @param drot Rotation increment matrix.
 * @param celent Characteristic element length, which is a typical length of a line across an element for a first-order element; it is half of the same typical length for a second-order element.
 * @param dfgrd0 Array containing the deformation gradient at the beginning of the increment.
 * @param dfgrd1 Array containing the deformation gradient at the end of the increment.
 * @param noel Element number.
 * @param npt Integration point number.
 * @param layer Layer number.
 * @param kspt Section point number within the current layer.
 * @param kstep Procedure type key; nlgeom boolean; linear perturbation boolean.
 * @param kinc Increment number.
 * */
extern "C" void umat(
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

#endif