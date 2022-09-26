#include <iostream>
#include <Eigen/Eigen>
#include "umat.cpp"

int main()
{      
    int ndi = 3;
    int nshr = 3;
    int ntens = ndi*nshr;
    double stress[ntens] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
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

    // Test Eigen library.
    int n = 10;
    Eigen::VectorXi v(n);
    for(int i = 0; i < n;  i++) {
        v[i] = i;
    }
    std::cout << "vector: " << v.transpose() << "\n";

    // Test mapping from C arrays to Eigen structures.
    int array[9];
    for(int i = 0; i < 9; ++i) {
        array[i] = i;
    }
    
    // Map to Eigen structure.
    typedef Eigen::Matrix<int, 3, 3> Matrix3i;
    Eigen::Map<Matrix3i> array_eigen(array);
    std::cout << array_eigen << "\n";
    
    // Modify Eigen atructure.
    array_eigen(0,1) = 100;
    std::cout << array_eigen << "\n";

    // Check this is reflected in the C array.
    for(int i = 0; i < 9; ++i) std::cout << array[i] << " ";
    std::cout << "\n";
}