
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "MCC.hpp"
#include "C2MC.hpp"

// MCC parameters: M=0.92, nu=0.2, N=1.195, lambda*=0.08, kappa*=0.02, p_c=10.
static Parameters mcc_params {{0.92, 0.2, 1.195, 0.08, 0.02}};

TEST_CASE("MCC::compute_f inside yield surface") {
    // p' = 5 (inside yield surface since p' < p_c), q = 0.
    // f = q^2 + M^2 * p' * (p' - p_c) = 0 + 0.92^2 * 5 * (5 - 10) = 0.8464 * (-25) = -21.16.
    State state {{10.0}};
    MCC model(mcc_params, state);

    Cauchy sigma_prime = Cauchy::Zero();
    sigma_prime(0,0) = 5.0;
    sigma_prime(1,1) = 5.0;
    sigma_prime(2,2) = 5.0;

    double f = model.mf.compute_f(sigma_prime, state);
    REQUIRE(f < 0.0);
}

TEST_CASE("MCC::compute_f on yield surface") {
    // At the apex of the yield surface: p' = p_c/2 = 5, q = M * p_c / 2 = 0.92 * 10 / 2 = 4.6.
    // f = q^2 + M^2 * p' * (p' - p_c) = 4.6^2 + 0.8464 * 5 * (5 - 10) = 21.16 - 21.16 = 0.
    State state {{10.0}};
    MCC model(mcc_params, state);

    double M = 0.92;
    double p_c = 10.0;
    double p_prime = p_c / 2.0;
    double q_target = M * p_c / 2.0;

    // Construct a stress state with the desired p' and q.
    // Use triaxial conditions: sigma_1 = p' + 2q/3, sigma_2 = sigma_3 = p' - q/3.
    Cauchy sigma_prime_tensor = Cauchy::Zero();
    sigma_prime_tensor(0,0) = p_prime + 2.0 * q_target / 3.0;
    sigma_prime_tensor(1,1) = p_prime - q_target / 3.0;
    sigma_prime_tensor(2,2) = p_prime - q_target / 3.0;

    double f = model.mf.compute_f(sigma_prime_tensor, state);
    REQUIRE_THAT(f, Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("MCC::solve elastic increment") {
    // Small isotropic strain inside the yield surface should not change state variables.
    State state {{10.0}};
    MCC model(mcc_params, state);

    // Isotropic stress state well inside the yield surface.
    Voigt sigma_tilde;
    sigma_tilde << 5.0, 5.0, 5.0, 0.0, 0.0, 0.0;
    model.set_sigma_prime_tilde(sigma_tilde);

    // Very small isotropic strain increment (stays elastic).
    Voigt dstran;
    dstran << 0.0001, 0.0001, 0.0001, 0.0, 0.0, 0.0;
    model.set_Delta_epsilon_tilde(dstran);

    model.solve();

    // State variables should remain unchanged for a purely elastic step.
    State result_state = model.get_state_variables();
    REQUIRE_THAT(result_state[0], Catch::Matchers::WithinRel(10.0, 0.001));
}

TEST_CASE("MCC::solve elastoplastic increment") {
    // Triaxial compression crossing the yield surface.
    State state {{10.0}};
    MCC model(mcc_params, state);

    // Start at isotropic stress equal to p_c (on the yield surface at q=0).
    Voigt sigma_tilde;
    sigma_tilde << 10.0, 10.0, 10.0, 0.0, 0.0, 0.0;
    model.set_sigma_prime_tilde(sigma_tilde);

    // Apply compressive strain in one direction to cross the yield surface.
    Voigt dstran;
    dstran << 0.005, -0.0025, -0.0025, 0.0, 0.0, 0.0;
    model.set_Delta_epsilon_tilde(dstran);

    model.solve();

    // After elastoplastic solve, f should be <= FTOL.
    Settings settings = model.get_settings();
    State result_state = model.get_state_variables();
    Voigt result_sigma = model.get_sigma_prime_tilde();
    Cauchy result_cauchy = to_cauchy(result_sigma);
    double f = model.mf.compute_f(result_cauchy, result_state);

    REQUIRE(f <= settings.FTOL);

    // State variables should have changed (p_c should have increased).
    REQUIRE(result_state[0] != 10.0);
}

TEST_CASE("MCC::solve Jacobian non-zero after elastoplastic solve") {
    // After an elastoplastic solve, the Jacobian should not be all zeros.
    State state {{10.0}};
    MCC model(mcc_params, state);

    Voigt sigma_tilde;
    sigma_tilde << 10.0, 10.0, 10.0, 0.0, 0.0, 0.0;
    model.set_sigma_prime_tilde(sigma_tilde);

    Voigt dstran;
    dstran << 0.005, -0.0025, -0.0025, 0.0, 0.0, 0.0;
    model.set_Delta_epsilon_tilde(dstran);

    model.solve();

    Jacobian J = model.get_jacobian();
    double norm = J.norm();
    REQUIRE(norm > 0.0);
}
