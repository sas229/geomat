
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "MCC.hpp"
#include <cmath>

// Use MCC model to provide test environment for Model base class methods.
static Parameters model_parameters {{0.92, 0.2, 1.195, 0.08, 0.02}};
static State model_state{{10.0}};

TEST_CASE("Model::compute_p") {
    MCC model(model_parameters, model_state);
    Cauchy sigma_prime = Cauchy::Zero();
    sigma_prime(0,0) = 100.0;
    sigma_prime(1,1) = 200.0;
    sigma_prime(2,2) = 300.0;
    double p_prime = model.compute_p_prime(sigma_prime);
    REQUIRE_THAT(p_prime, Catch::Matchers::WithinRel(200.0, 0.001));
}

TEST_CASE("Model::compute_q") {
    MCC model(model_parameters, model_state);
    Cauchy sigma_prime = Cauchy::Zero();
    sigma_prime(0,0) = 49.0;
    sigma_prime(1,1) = 50.0;
    sigma_prime(2,2) = 50.0;
    double q = model.compute_q(sigma_prime);
    REQUIRE_THAT(q, Catch::Matchers::WithinRel(1.0, 0.001));
}

TEST_CASE("Model::compute_stress_invariants") {
    MCC model(model_parameters, model_state);
    // Set up a known diagonal stress state (no shear) for simpler verification.
    Cauchy sigma = Cauchy::Zero();
    sigma(0,0) = 100.0;
    sigma(1,1) = 200.0;
    sigma(2,2) = 300.0;

    double I_1, I_2, I_3, J_1, J_2, J_3;
    model.compute_stress_invariants(sigma, I_1, I_2, I_3, J_1, J_2, J_3);

    // I_1 = trace(sigma) = 100 + 200 + 300 = 600.
    REQUIRE_THAT(I_1, Catch::Matchers::WithinRel(600.0, 0.001));

    // I_2 = 0.5 * (trace(sigma)^2 - trace(sigma^2)).
    // trace(sigma^2) = 100^2 + 200^2 + 300^2 = 140000.
    // I_2 = 0.5 * (360000 - 140000) = 110000.
    REQUIRE_THAT(I_2, Catch::Matchers::WithinRel(110000.0, 0.001));

    // I_3 = det(sigma) = 100 * 200 * 300 = 6000000.
    REQUIRE_THAT(I_3, Catch::Matchers::WithinRel(6000000.0, 0.001));

    // J_1 should always be zero (trace of deviatoric tensor).
    REQUIRE_THAT(J_1, Catch::Matchers::WithinAbs(0.0, 1e-6));

    // J_2 >= 0.
    REQUIRE(J_2 >= 0.0);

    // Verify J_2 and J_3 consistency against manually computed deviatoric tensor.
    double p = I_1 / 3.0;
    Cauchy s = sigma - p * Cauchy::Identity();
    double J_2_check = 0.5 * (s * s).trace();
    REQUIRE_THAT(J_2, Catch::Matchers::WithinRel(J_2_check, 0.001));

    double J_3_check = s.determinant();
    REQUIRE_THAT(J_3, Catch::Matchers::WithinRel(J_3_check, 0.001));
}

TEST_CASE("Model::compute_principal_stresses") {
    MCC model(model_parameters, model_state);
    // Diagonal stress tensor: eigenvalues should match diagonal entries.
    Cauchy sigma_prime = Cauchy::Zero();
    sigma_prime(0,0) = 300.0;
    sigma_prime(1,1) = 200.0;
    sigma_prime(2,2) = 100.0;

    double sigma_1, sigma_2, sigma_3;
    Cauchy R, S;
    model.compute_principal_stresses(sigma_prime, sigma_1, sigma_2, sigma_3, R, S);

    // sigma_1 is the maximum eigenvalue, sigma_3 is the minimum.
    REQUIRE_THAT(sigma_1, Catch::Matchers::WithinRel(300.0, 0.001));
    REQUIRE_THAT(sigma_2, Catch::Matchers::WithinRel(200.0, 0.001));
    REQUIRE_THAT(sigma_3, Catch::Matchers::WithinRel(100.0, 0.001));
}

TEST_CASE("Model::compute_lode") {
    MCC model(model_parameters, model_state);
    // Triaxial compression state: sigma_1 > sigma_2 = sigma_3.
    // With sigma_1 = 200, sigma_2 = sigma_3 = 100, the deviatoric stress has
    // s_1 > 0, s_2 = s_3 < 0, giving J_3 > 0 and theta_s_bar = -pi/6.
    Cauchy sigma = Cauchy::Zero();
    sigma(0,0) = 200.0;
    sigma(1,1) = 100.0;
    sigma(2,2) = 100.0;

    double I_1, I_2, I_3, J_1, J_2, J_3;
    model.compute_stress_invariants(sigma, I_1, I_2, I_3, J_1, J_2, J_3);

    double theta_c, theta_s, theta_s_bar;
    model.compute_lode(J_2, J_3, theta_c, theta_s, theta_s_bar);

    // Verify theta_s_bar is within the valid range [-pi/6, pi/6].
    REQUIRE(theta_s_bar >= -pi / 6.0 - 1e-6);
    REQUIRE(theta_s_bar <= pi / 6.0 + 1e-6);

    // For this stress state, theta_s_bar should be near -pi/6.
    REQUIRE_THAT(theta_s_bar, Catch::Matchers::WithinAbs(-pi / 6.0, 1e-6));
}

TEST_CASE("Model::compute_s") {
    MCC model(model_parameters, model_state);
    Cauchy sigma = Cauchy::Zero();
    sigma(0,0) = 100.0;
    sigma(1,1) = 200.0;
    sigma(2,2) = 300.0;

    double p = model.compute_p(sigma);
    Cauchy s = model.compute_s(sigma, p);

    // The deviatoric tensor should have zero trace.
    double trace_s = s.trace();
    REQUIRE_THAT(trace_s, Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("Model::compute_mises_stress") {
    MCC model(model_parameters, model_state);
    // For a known J_2, the Mises stress should be sqrt(3 * J_2).
    double J_2 = 100.0;
    double mises = model.compute_mises_stress(J_2);
    double expected = std::sqrt(3.0 * J_2);
    REQUIRE_THAT(mises, Catch::Matchers::WithinRel(expected, 0.001));
}

TEST_CASE("Model::compute_max_shear") {
    MCC model(model_parameters, model_state);
    // max_shear = (sigma_1 - sigma_3) / 2.
    double sigma_1 = 300.0;
    double sigma_2 = 200.0;
    double sigma_3 = 100.0;
    double max_shear = model.compute_max_shear(sigma_1, sigma_2, sigma_3);
    double expected = (sigma_1 - sigma_3) / 2.0;
    REQUIRE_THAT(max_shear, Catch::Matchers::WithinRel(expected, 0.001));
}

TEST_CASE("Model::compute_Delta_epsilon_vol") {
    MCC model(model_parameters, model_state);
    // Volumetric strain increment = trace of strain increment tensor.
    Cauchy Delta_epsilon = Cauchy::Zero();
    Delta_epsilon(0,0) = 0.001;
    Delta_epsilon(1,1) = 0.002;
    Delta_epsilon(2,2) = 0.003;

    double vol = model.compute_Delta_epsilon_vol(Delta_epsilon);
    double expected = 0.001 + 0.002 + 0.003;
    REQUIRE_THAT(vol, Catch::Matchers::WithinRel(expected, 0.001));
}
