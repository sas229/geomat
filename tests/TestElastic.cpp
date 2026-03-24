
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "LinearElastic.hpp"

TEST_CASE("LinearElastic::compute_D_e") {
    // Initialise LinearElastic with K=100, G=50.
    Parameters params(2);
    params << 100.0, 50.0;
    State state = State::Zero(0);
    LinearElastic model(params, state);

    // Compute elastic stiffness matrix.
    Cauchy sigma_prime = Cauchy::Zero();
    Constitutive D_e = model.compute_isotropic_linear_elastic_matrix(100.0, 50.0);

    // D_e(0,0) = K + 4/3*G = 100 + 200/3 = 166.6667.
    double expected_00 = 100.0 + 4.0 / 3.0 * 50.0;
    REQUIRE_THAT(D_e(0,0), Catch::Matchers::WithinRel(expected_00, 0.001));

    // D_e(0,1) = K - 2/3*G = 100 - 100/3 = 66.6667.
    double expected_01 = 100.0 - 2.0 / 3.0 * 50.0;
    REQUIRE_THAT(D_e(0,1), Catch::Matchers::WithinRel(expected_01, 0.001));

    // D_e(3,3) = G = 50.
    REQUIRE_THAT(D_e(3,3), Catch::Matchers::WithinRel(50.0, 0.001));
}

TEST_CASE("LinearElastic::solve stress increment") {
    // Initialise LinearElastic with K=100, G=50.
    Parameters params(2);
    params << 100.0, 50.0;
    State state = State::Zero(0);
    LinearElastic model(params, state);

    // Set initial stress to zero.
    Voigt sigma_tilde = Voigt::Zero();
    model.set_sigma_prime_tilde(sigma_tilde);

    // Apply a strain increment.
    Voigt dstran;
    dstran << 0.001, 0.0, 0.0, 0.0, 0.0, 0.0;
    model.set_Delta_epsilon_tilde(dstran);

    // Solve.
    model.solve();

    // Expected stress = D_e * dstran.
    Constitutive D_e = model.compute_isotropic_linear_elastic_matrix(100.0, 50.0);
    Voigt expected_stress = D_e * dstran;

    Voigt result_stress = model.get_sigma_prime_tilde();
    for (int i = 0; i < 6; i++) {
        REQUIRE_THAT(result_stress[i], Catch::Matchers::WithinRel(expected_stress[i], 0.001));
    }
}

TEST_CASE("LinearElastic::solve Jacobian equals D_e") {
    // After solving, the Jacobian should equal the elastic matrix.
    Parameters params(2);
    params << 100.0, 50.0;
    State state = State::Zero(0);
    LinearElastic model(params, state);

    Voigt sigma_tilde = Voigt::Zero();
    model.set_sigma_prime_tilde(sigma_tilde);

    Voigt dstran;
    dstran << 0.001, 0.0, 0.0, 0.0, 0.0, 0.0;
    model.set_Delta_epsilon_tilde(dstran);
    model.solve();

    Jacobian J = model.get_jacobian();
    Constitutive D_e = model.compute_isotropic_linear_elastic_matrix(100.0, 50.0);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            REQUIRE_THAT(J(i,j), Catch::Matchers::WithinRel(D_e(i,j), 0.001));
        }
    }
}
