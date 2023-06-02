#include <iostream>
#include <typeinfo>

#include "umat.hpp"

int main() {      
    // int ndi = 3;
    // int nshr = 3;
    // int ntens = 6;
    // double stress[ntens] = {-1000, -1000, -1000, 0, 0, 0};
    // int nstatv = 0;
    // double statev[nstatv] = {};
    // double ddsdde[ntens*ntens] = {0};
    // double sse, spd, scd, rpl;
    // double ddsddt[ntens], drplde[ntens];
    // double drpldt;
    // double stran[ntens] = {0};
    // double dstran[ntens] = {-0.01, -0.005, -0.005, 0.0, 0.0, 0.0};
    // double time[2];
    // double dtime, temp, dtemp;
    // double predef, dpred;
    // char cmname[] = "LinearElastic";
    // int nprops = 2;
    // double props[nprops] = {50000, 25000};
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

    // int ndi = 3;
    // int nshr = 3;
    // int ntens = 6;
    // double stress[ntens] = {-10, -10, -10, 0, 0, 0};
    // int nstatv = 2;
    // double statev[nstatv] = {1.7477796692480023, 10};
    // double ddsdde[ntens*ntens] = {0};
    // double sse, spd, scd, rpl;
    // double ddsddt[ntens], drplde[ntens];
    // double drpldt;
    // double stran[ntens] = {0};
    // double dstran[ntens] = {-0.0005, 0.00025, 0.00025, 0.0, 0.0, 0.0};
    // double time[2];
    // double dtime, temp, dtemp;
    // double predef, dpred;
    // char cmname[] = "MCC";
    // int nprops = 5;
    // double props[nprops] = {0.92, 0.2, 1.195, 0.08, 0.02};
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

    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    double stress[ntens] = {-10, -10, -10, 0, 0, 0};
    int nstatv = 3;
    double statev[nstatv] = {1.7477796692480023, 10, 3.0};
    double ddsdde[ntens*ntens] = {0};
    double sse, spd, scd, rpl;
    double ddsddt[ntens], drplde[ntens];
    double drpldt;
    double stran[ntens] = {0};
    double dstran[ntens] = {-0.0005, 0.00025, 0.00025, 0.0, 0.0, 0.0};
    double time[2];
    double dtime, temp, dtemp;
    double predef, dpred;
    char cmname[] = "SMCC";
    int nprops = 8;
    double props[nprops] = {0.92, 0.2, 1.195, 0.08, 0.02, 3.0, 0.2, 0.5};
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
    for (int i=0; i < 100; ++i) {
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