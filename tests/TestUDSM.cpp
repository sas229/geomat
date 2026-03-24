
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "udsm_cpp_interface.hpp"
#include "udsm_metadata.hpp"

// Helper to initialise common UDSM call arguments.
struct UDSMArgs {
    int IDTask = 0, iMod = 1, IsUndr = 0, iStep = 0, iTer = 0, iEl = 1, Int_val = 1;
    double X = 0, Y = 0, Z = 0, Time0 = 0, dTime = 0;
    double Props[50] = {0};
    double Sig0[6] = {0};
    double Swp0 = 0;
    double StVar0[10] = {0};
    double dEps[6] = {0};
    double D[36] = {0};
    double Bulk_W = 0;
    double Sig[6] = {0};
    double Swp = 0;
    double StVar[10] = {0};
    int ipl = 0, nStat = 0, NonSym = 0, iStrsDep = 0, iTimeDep = 0, iTang = 0;
    char iPrjDir[256] = {0};
    int iPrjLen = 0;
    int iAbort = 0;

    void call() {
        udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int_val,
            &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
            dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
            &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);
    }
};

TEST_CASE("UDSM IDTask 4 returns correct nStat") {
    // LinearElastic (iMod=1): 0 state variables.
    {
        UDSMArgs args;
        args.IDTask = 4;
        args.iMod = 1;
        args.call();
        REQUIRE(args.nStat == 0);
    }

    // MCC (iMod=4): 1 state variable.
    {
        UDSMArgs args;
        args.IDTask = 4;
        args.iMod = 4;
        args.call();
        REQUIRE(args.nStat == 1);
    }

    // SMCC (iMod=5): 2 state variables.
    {
        UDSMArgs args;
        args.IDTask = 4;
        args.iMod = 5;
        args.call();
        REQUIRE(args.nStat == 2);
    }
}

TEST_CASE("UDSM IDTask 5 returns correct matrix attributes") {
    UDSMArgs args;
    args.IDTask = 5;
    args.call();

    REQUIRE(args.NonSym == 0);
    REQUIRE(args.iStrsDep == 1);
    REQUIRE(args.iTimeDep == 0);
    REQUIRE(args.iTang == 0);
}

TEST_CASE("UDSM IDTask 6 elastic stiffness for LinearElastic") {
    // LinearElastic (iMod=1) with K=100, G=50.
    UDSMArgs args;
    args.IDTask = 6;
    args.iMod = 1;
    args.Props[0] = 100.0;
    args.Props[1] = 50.0;

    // Set a non-zero stress state (needed for some models, harmless for LinearElastic).
    args.Sig0[0] = -10.0;
    args.Sig0[1] = -10.0;
    args.Sig0[2] = -10.0;

    args.call();

    // D is stored in column-major (Fortran) order: D[j*6 + i] = D_e(i,j).
    double K = 100.0;
    double G = 50.0;
    double expected_D00 = K + 4.0 / 3.0 * G;
    REQUIRE_THAT(args.D[0 * 6 + 0], Catch::Matchers::WithinRel(expected_D00, 0.001));

    // D_e(0,1) = K - 2/3*G.
    double expected_D01 = K - 2.0 / 3.0 * G;
    REQUIRE_THAT(args.D[1 * 6 + 0], Catch::Matchers::WithinRel(expected_D01, 0.001));

    // D_e(3,3) = G.
    REQUIRE_THAT(args.D[3 * 6 + 3], Catch::Matchers::WithinRel(G, 0.001));
}

TEST_CASE("UDSM GetModelCount returns 5") {
    int count = 0;
    GetModelCount(&count);
    REQUIRE(count == 5);
}

TEST_CASE("UDSM GetParamCount for each model") {
    // Expected parameter counts: LinearElastic=2, C2MC=7, EMC=8, MCC=5, SMCC=8.
    int expected_counts[5] = {2, 7, 8, 5, 8};
    for (int i = 1; i <= 5; i++) {
        int count = 0;
        GetParamCount(&i, &count);
        REQUIRE(count == expected_counts[i - 1]);
    }
}
