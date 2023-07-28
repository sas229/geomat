
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "MCC.hpp"
#include <iostream>

// Use MCC model to provide test environment.
Parameters parameters {{0.92, 0.2, 1.195, 0.08, 0.02}};
State state{{10.0}}; 
Cauchy sigma_prime = Cauchy::Zero(); 
MCC model(parameters, state);

TEST_CASE("Model::compute_p") {
    sigma_prime(0,0) = 100.0;
    sigma_prime(1,1) = 200.0;
    sigma_prime(2,2) = 300.0;
    double p_prime = model.compute_p_prime(sigma_prime);
    REQUIRE_THAT(p_prime, Catch::Matchers::WithinRel(200.0, 0.001));
}

TEST_CASE("Model::compute_q") {
    sigma_prime(0,0) = 49.0;
    sigma_prime(1,1) = 50.0;
    sigma_prime(2,2) = 50.0;
    double q = model.compute_q(sigma_prime);
    REQUIRE_THAT(q, Catch::Matchers::WithinRel(1.0, 0.001));
}