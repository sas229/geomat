
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "umat_cpp_interface.hpp"

TEST_CASE("umat_cpp_interface LinearElastic stress increment") {
    // Set up LinearElastic model with K=100, G=50.
    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    int nstatv = 0;
    int nprops = 2;
    double props[2] = {100.0, 50.0};

    // Abaqus uses tension-positive convention; internal model uses compression-positive.
    // Initial stress is zero.
    double stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double statev_dummy = 0.0;
    double ddsdde[36] = {0};
    double sse = 0, spd = 0, scd = 0, rpl = 0;
    double ddsddt[6] = {0};
    double drplde[6] = {0};
    double drpldt = 0;
    double stran[6] = {0};
    // Tension-positive strain increment in Abaqus convention.
    double dstran[6] = {0.001, 0.0, 0.0, 0.0, 0.0, 0.0};
    double time_arr[2] = {0, 0};
    double dtime = 1.0;
    double temp_val = 0, dtemp = 0;
    double predef = 0, dpred = 0;
    char cmname[] = "LinearElastic";
    double coords[3] = {0};
    double drot = 0;
    double pnewdt = 0, celent = 0;
    double dfgrd0 = 0, dfgrd1 = 0;
    int noel = 1, npt = 1, layer = 0, kspt = 0, kstep = 0, kinc = 1;

    umat_cpp_interface(
        stress, &statev_dummy, ddsdde, &sse, &spd, &scd, &rpl,
        ddsddt, drplde, &drpldt, stran, dstran, time_arr, &dtime,
        &temp_val, &dtemp, &predef, &dpred, cmname, &ndi, &nshr, &ntens,
        &nstatv, props, &nprops, coords, &drot, &pnewdt, &celent,
        &dfgrd0, &dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc
    );

    // Expected: D_e * (-dstran) negated back = D_e * dstran in tension-positive.
    // D_e(0,0) = K + 4/3*G = 100 + 200/3.
    // stress[0] should be -(-(K+4/3*G)*0.001) = (K+4/3*G)*0.001.
    double K = 100.0;
    double G = 50.0;
    double expected_s0 = (K + 4.0/3.0*G) * 0.001;
    double expected_s1 = (K - 2.0/3.0*G) * 0.001;

    // Abaqus convention: the result stress should match D_e * dstran (tension-positive).
    REQUIRE_THAT(stress[0], Catch::Matchers::WithinRel(expected_s0, 0.001));
    REQUIRE_THAT(stress[1], Catch::Matchers::WithinRel(expected_s1, 0.001));
    REQUIRE_THAT(stress[2], Catch::Matchers::WithinRel(expected_s1, 0.001));
    REQUIRE_THAT(stress[3], Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("umat_cpp_interface MCC state update") {
    // Run one increment with MCC and verify state variables are updated.
    int ndi = 3;
    int nshr = 3;
    int ntens = 6;
    int nstatv = 1;
    int nprops = 5;
    double props[5] = {0.92, 0.2, 1.195, 0.08, 0.02};

    // Abaqus tension-positive: compression is negative.
    double stress[6] = {-10.0, -10.0, -10.0, 0.0, 0.0, 0.0};
    double statev[1] = {10.0};
    double ddsdde[36] = {0};
    double sse = 0, spd = 0, scd = 0, rpl = 0;
    double ddsddt[6] = {0};
    double drplde[6] = {0};
    double drpldt = 0;
    double stran[6] = {0};
    // Compression in Abaqus = negative strain.
    double dstran[6] = {-0.005, 0.0025, 0.0025, 0.0, 0.0, 0.0};
    double time_arr[2] = {0, 0};
    double dtime = 1.0;
    double temp_val = 0, dtemp = 0;
    double predef = 0, dpred = 0;
    char cmname[] = "MCC";
    double coords[3] = {0};
    double drot = 0;
    double pnewdt = 0, celent = 0;
    double dfgrd0 = 0, dfgrd1 = 0;
    int noel = 1, npt = 1, layer = 0, kspt = 0, kstep = 0, kinc = 1;

    umat_cpp_interface(
        stress, statev, ddsdde, &sse, &spd, &scd, &rpl,
        ddsddt, drplde, &drpldt, stran, dstran, time_arr, &dtime,
        &temp_val, &dtemp, &predef, &dpred, cmname, &ndi, &nshr, &ntens,
        &nstatv, props, &nprops, coords, &drot, &pnewdt, &celent,
        &dfgrd0, &dfgrd1, &noel, &npt, &layer, &kspt, &kstep, &kinc
    );

    // State variable (p_c) should have been updated.
    REQUIRE(statev[0] != 10.0);
}
