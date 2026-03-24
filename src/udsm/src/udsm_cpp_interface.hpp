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
