#include <iostream>
#include <cstring>

#include "udsm_cpp_interface.hpp"
#include "udsm_metadata.hpp"

int main() {
    // Test with MCC model (iMod=4 in the registry: LinearElastic=1, C2MC=2, EMC=3, MCC=4, SMCC=5).
    int iMod = 4;

    // ---------------------------
    // IDTask 4: Query state variable count.
    // ---------------------------
    int IDTask = 4;
    int nStat = 0;
    int IsUndr = 0, iStep = 0, iTer = 0, iEl = 1, Int = 1;
    double X = 0, Y = 0, Z = 0, Time0 = 0, dTime = 0;
    double Props[50] = {0};
    double Sig0[20] = {0};
    double Swp0 = 0;
    double StVar0[10] = {0};
    double dEps[12] = {0};
    double D[36] = {0};
    double Bulk_W = 0;
    double Sig[6] = {0};
    double Swp = 0;
    double StVar[10] = {0};
    int ipl = 0, NonSym = 0, iStrsDep = 0, iTimeDep = 0, iTang = 0;
    char iPrjDir[256] = {0};
    int iPrjLen = 0;
    int iAbort = 0;

    udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
        &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
        dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
        &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

    std::cout << "IDTask 4 - nStat: " << nStat << " (expected 1)" << std::endl;

    // ---------------------------
    // IDTask 5: Query matrix attributes.
    // ---------------------------
    IDTask = 5;
    udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
        &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
        dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
        &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

    std::cout << "IDTask 5 - NonSym: " << NonSym << ", iStrsDep: " << iStrsDep
              << ", iTimeDep: " << iTimeDep << ", iTang: " << iTang << std::endl;

    // ---------------------------
    // Set up MCC model parameters: M=0.92, nu=0.2, N=1.195, lambda*=0.08, kappa*=0.02
    // ---------------------------
    Props[0] = 0.92;
    Props[1] = 0.2;
    Props[2] = 1.195;
    Props[3] = 0.08;
    Props[4] = 0.02;

    // Initial stress state (compression-negative for Plaxis): p' = 7.5 kPa isotropic.
    Sig0[0] = -7.5;  // sigma_xx
    Sig0[1] = -7.5;  // sigma_yy
    Sig0[2] = -7.5;  // sigma_zz
    Sig0[3] = 0.0;   // sigma_xy
    Sig0[4] = 0.0;   // sigma_yz
    Sig0[5] = 0.0;   // sigma_zx

    // Initial state variable: preconsolidation pressure p_c = 10 kPa.
    StVar0[0] = 10.0;

    // ---------------------------
    // IDTask 1: Initialise state variables.
    // ---------------------------
    IDTask = 1;
    udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
        &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
        dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
        &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

    // Override StVar0 with our initial p_c (as a user would in Plaxis).
    StVar0[0] = 10.0;
    std::cout << "IDTask 1 - State variables initialised. StVar0[0] = " << StVar0[0] << std::endl;

    // ---------------------------
    // IDTask 6: Compute elastic stiffness matrix.
    // ---------------------------
    IDTask = 6;
    std::memset(D, 0, sizeof(D));
    udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
        &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
        dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
        &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

    std::cout << "IDTask 6 - Elastic stiffness matrix D:" << std::endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            std::cout << D[j * 6 + i] << "\t";
        }
        std::cout << std::endl;
    }

    // ---------------------------
    // IDTask 2: Calculate constitutive stresses (run increments).
    // ---------------------------
    dEps[0] = -0.0005;  // Δε_xx (compression)
    dEps[1] = 0.00025;  // Δε_yy
    dEps[2] = 0.00025;  // Δε_zz

    std::cout << "\nIDTask 2 - Running 100 strain increments:" << std::endl;
    std::cout << "Sig_xx\tSig_yy\tSig_zz\tSig_xy\tSig_yz\tSig_zx" << std::endl;

    for (int inc = 0; inc < 100; inc++) {
        IDTask = 2;
        udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
            &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
            dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
            &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

        // Update Sig0 and StVar0 for next increment.
        for (int i = 0; i < 6; i++) Sig0[i] = Sig[i];
        for (int i = 0; i < nStat; i++) StVar0[i] = StVar[i];

        std::cout << Sig[0] << "\t" << Sig[1] << "\t" << Sig[2] << "\t"
                  << Sig[3] << "\t" << Sig[4] << "\t" << Sig[5] << std::endl;
    }

    // ---------------------------
    // IDTask 3: Compute effective stiffness matrix.
    // ---------------------------
    IDTask = 3;
    std::memset(D, 0, sizeof(D));
    udsm_cpp_interface(&IDTask, &iMod, &IsUndr, &iStep, &iTer, &iEl, &Int,
        &X, &Y, &Z, &Time0, &dTime, Props, Sig0, &Swp0, StVar0,
        dEps, D, &Bulk_W, Sig, &Swp, StVar, &ipl, &nStat, &NonSym,
        &iStrsDep, &iTimeDep, &iTang, iPrjDir, &iPrjLen, &iAbort);

    std::cout << "\nIDTask 3 - Effective stiffness matrix D:" << std::endl;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            std::cout << D[j * 6 + i] << "\t";
        }
        std::cout << std::endl;
    }

    // ---------------------------
    // Test metadata functions.
    // ---------------------------
    std::cout << "\n--- Metadata Functions ---" << std::endl;

    int modelCount = 0;
    GetModelCount(&modelCount);
    std::cout << "GetModelCount: " << modelCount << std::endl;

    for (int m = 1; m <= modelCount; m++) {
        char name[256] = {0};
        GetModelName(&m, name);
        int nameLen = static_cast<unsigned char>(name[0]);
        std::cout << "Model " << m << ": " << std::string(name + 1, nameLen);

        int paramCount = 0;
        GetParamCount(&m, &paramCount);
        std::cout << " (" << paramCount << " params";

        int stateCount = 0;
        GetStateVarCount(&m, &stateCount);
        std::cout << ", " << stateCount << " state vars)" << std::endl;

        for (int p = 1; p <= paramCount; p++) {
            char pname[256] = {0};
            char punit[256] = {0};
            GetParamName(&m, &p, pname);
            GetParamUnit(&m, &p, punit);
            int pnameLen = static_cast<unsigned char>(pname[0]);
            int punitLen = static_cast<unsigned char>(punit[0]);
            std::cout << "  Param " << p << ": " << std::string(pname + 1, pnameLen)
                      << " [" << std::string(punit + 1, punitLen) << "]" << std::endl;
        }

        for (int s = 1; s <= stateCount; s++) {
            char sname[256] = {0};
            char sunit[256] = {0};
            GetStateVarName(&m, &s, sname);
            GetStateVarUnit(&m, &s, sunit);
            int snameLen = static_cast<unsigned char>(sname[0]);
            int sunitLen = static_cast<unsigned char>(sunit[0]);
            std::cout << "  State " << s << ": " << std::string(sname + 1, snameLen)
                      << " [" << std::string(sunit + 1, sunitLen) << "]" << std::endl;
        }
    }

    std::cout << "\nAll UDSM integration tests completed." << std::endl;
    return 0;
}
