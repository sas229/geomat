
#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Tensor.hpp"

TEST_CASE("Tensor::to_voigt and to_cauchy round trip") {
    // Create a Cauchy tensor with known values.
    Cauchy original = Cauchy::Zero();
    original(0,0) = 100.0;
    original(1,1) = 200.0;
    original(2,2) = 300.0;
    original(0,1) = 10.0;
    original(1,0) = 10.0;
    original(0,2) = 20.0;
    original(2,0) = 20.0;
    original(1,2) = 30.0;
    original(2,1) = 30.0;

    // Round trip: Cauchy -> Voigt -> Cauchy.
    Voigt voigt = to_voigt(original);
    Cauchy recovered = to_cauchy(voigt);

    // Verify all components match.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            REQUIRE_THAT(recovered(i,j), Catch::Matchers::WithinRel(original(i,j), 0.001));
        }
    }
}

TEST_CASE("Tensor::dev has zero trace") {
    Cauchy tensor = Cauchy::Zero();
    tensor(0,0) = 100.0;
    tensor(1,1) = 200.0;
    tensor(2,2) = 300.0;
    tensor(0,1) = 10.0;
    tensor(1,0) = 10.0;

    Cauchy deviatoric = dev(tensor);
    double trace_dev = tr(deviatoric);

    REQUIRE_THAT(trace_dev, Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("Tensor::double_dot_product of identity tensor") {
    // The double dot product (Frobenius inner product) of the identity with itself
    // should equal 3.0 (sum of squares of all entries: 1^2 + 1^2 + 1^2 = 3).
    Cauchy eye = Cauchy::Identity();
    double result = double_dot_product(eye);

    REQUIRE_THAT(result, Catch::Matchers::WithinRel(3.0, 0.001));
}

TEST_CASE("Tensor::double_dot_product of two tensors") {
    Cauchy a = Cauchy::Identity();
    Cauchy b = Cauchy::Identity();
    double result = double_dot_product(a, b);

    REQUIRE_THAT(result, Catch::Matchers::WithinRel(3.0, 0.001));
}

TEST_CASE("Tensor::tr equals sum of diagonal") {
    Cauchy tensor = Cauchy::Zero();
    tensor(0,0) = 10.0;
    tensor(1,1) = 20.0;
    tensor(2,2) = 30.0;
    tensor(0,1) = 5.0;
    tensor(1,0) = 5.0;

    double trace = tr(tensor);
    double expected = 10.0 + 20.0 + 30.0;

    REQUIRE_THAT(trace, Catch::Matchers::WithinRel(expected, 0.001));
}
