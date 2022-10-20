#include <iostream>
#include <typeinfo>

#include "umat.hpp"

int main() {      
    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    double stress[ntens] = {-1000, -1000, -1000, 0, 0, 0};
    int nstatv = 0;
    double statev[nstatv] = {};
    double ddsdde[ndi*nshr] = {};
    double sse, spd, scd, rpl;
    double ddsddt[ntens], drplde[ntens];
    double drpldt;
    double stran[ntens] = {0};
    double dstran[ntens] = {-0.01, -0.005, -0.005, 0.0, 0.0, 0.0};
    double time[2];
    double dtime, temp, dtemp;
    double predef, dpred;
    char cmname[] = "LinearElastic";
    int nprops = 2;
    double props[nprops] = {50000, 25000};
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

    // int ndi = 3;
    // int nshr = 3;
    // int ntens = 6;
    // double stress[ntens] = {-1000, -1000, -1000, 0, 0, 0};
    // int nstatv = 2;
    // double statev[nstatv] = {2.0, 1200};
    // double ddsdde[ndi*nshr] = {3, 2, 1, 6, 5, 4, 9, 8, 7};
    // double sse, spd, scd, rpl;
    // double ddsddt[ntens], drplde[ntens];
    // double drpldt;
    // double stran[ntens] = {0};
    // double dstran[ntens] = {-0.01, -0.005, -0.005, 0.0, 0.0, 0.0};
    // double time[2];
    // double dtime, temp, dtemp;
    // double predef, dpred;
    // char cmname[] = "MCC";
    // int nprops = 5;
    // double props[nprops] = {0.9, 0.3, 1.10, 0.205, 0.044};
    // double coords[3] = {0.0, 0.0, 0.0};
    // double drot;
    // double pnewdt, celent;
    // double dfgrd0;
    // double dfgrd1;
    // int noel = 1;
    // int npt = 1;
    // int layer;
    // int kspt;
    // int kstep;
    // int kinc;

    // Try ten increments.
    for (int i=0; i < 1; ++i) {
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