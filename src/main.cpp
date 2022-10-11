#include <iostream>
#include <typeinfo>

#include "umat.hpp"

int main() {      
    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    double stress[ntens] = {1, 2, 3, 4, 5, 6};
    int nstatv = 2;
    double statev[nstatv] = {7, 8};
    double ddsdde[ndi*nshr] = {3, 2, 1, 6, 5, 4, 9, 8, 7};
    double sse, spd, scd, rpl;
    double ddsddt[ntens], drplde[ntens];
    double drpldt;
    double stran[ntens] = {0};
    double dstran[ntens] = {0.1, 0.1, 0.1, 0.0, 0.0, 0.0};
    double time[2];
    double dtime, temp, dtemp;
    double predef, dpred;
    char cmname[] = "MCC";
    int nprops = 2;
    double props[nprops];
    double coords[3] = {0.0, 0.0, 0.0};
    double drot;
    double pnewdt, celent;
    double dfgrd0;
    double dfgrd1;
    int noel = 1;
    int npt = 1;
    int layer;
    int kspt;
    int kstep;
    int kinc;

    // Try ten increments.
    for (int i=0; i < 10; ++i) {
        umat(
            stress,
            statev,
            ddsdde,
            &sse,
            &spd,
            &scd,
            &rpl,
            ddsddt,
            drplde,
            &drpldt,
            stran,
            dstran,
            time,
            &dtime,
            &temp,
            &dtemp,
            &predef,
            &dpred,
            cmname,
            &ndi,
            &nshr,
            &ntens,
            &nstatv,
            props,
            &nprops,
            coords,
            &drot,
            &pnewdt,
            &celent,
            &dfgrd0,
            &dfgrd1,
            &noel,
            &npt,
            &layer,
            &kspt,
            &kstep,
            &kinc
        );
    }
}