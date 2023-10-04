#include <aba_for_c.h>
#include <iostream>

extern "C" void umat(double *stress, double *statev, double   *ddsdde,
      double *sse, double *spd, double *scd, double *rpl, double *ddsddt,
      double *drplde, double *drpldt, double *stran, double *dstran,
      double *time, double *dtime, double *temp, double *dtemp, double *predef,
      double *dpred, char *cmname, int *ndi,
      int *nshr, int *ntens, int *nstatv, double *props,
      int *nprops, double *coords, double drot[][3], double *pnewdt,
      double *celent, double  dfgrd0[][3], double dfgrd1[][3], int *noel,
      int *npt,  int *layer,  int *kspt,  int *kstep,
      int *kinc)
{
    std::cout << "Inside the UMAT from C++..." << std::endl;
}