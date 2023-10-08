#include <iostream>
#include <typeinfo>

#include "umat_cpp_interface.hpp"

int main() {   
    // Test using the MCC model.
    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    double stress[ntens] = {-10, -10, -10, 0, 0, 0};
    int nstatv = 1;
    double statev[nstatv] = {10};
    double ddsdde[ntens*ntens] = {0};
    double sse, spd, scd, rpl;
    double ddsddt[ntens], drplde[ntens];
    double drpldt;
    double stran[ntens] = {0};
    double dstran[ntens] = {-0.0005, 0.00025, 0.00025, 0.0, 0.0, 0.0};
    double time[2];
    double dtime, temp, dtemp;
    double predef, dpred;
    char cmname[] = "MCC";
    int nprops = 5;
    double props[nprops] = {0.92, 0.2, 1.195, 0.08, 0.02};
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

    // // Test using C2mC model.   
    // int ndi = 3;
    // int nshr = 3;
    // int ntens = 6;
    // double stress[6] = {-50, -50, -50, 0, 0, 0};
    // int nstatv = 0;
    // double *statev; // Need to declare a pointer to a double because Windows cannot initialise an array of zero length!
    // double ddsdde[36] = {0};
    // double sse, spd, scd, rpl;
    // double ddsddt[6], drplde[6];
    // double drpldt;
    // double stran[6] = {0};
    // double dstran[6] = {-0.005, 0.0025, 0.0025, 0.0, 0.0, 0.0};
    // double time[2];
    // double dtime, temp, dtemp;
    // double predef, dpred;
    // char cmname[] = "C2MC";
    // int nprops = 7;
    // double props[7] = {500.0, 0.2, 0.0, 33.0, 5.0, 29.0, 1.1};
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

    // Print stress state.
    for (auto i=0; i<6; i++) {
        std::cout << stress[i] << " ";
    }
    std::cout << "\n";

    // Try some increments.
    for (int i=0; i < 100; ++i) {
        umat_cpp_interface(
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
        // Print stress state.
        for (auto i=0; i<6; i++) {
            std::cout << stress[i] << " ";
        }
        std::cout << "\n";
    }
}