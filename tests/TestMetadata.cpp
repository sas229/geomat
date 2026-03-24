
#include <catch2/catch_all.hpp>
#include "ModelLookup.hpp"
#include "Model.hpp"
#include "models.hpp"

TEST_CASE("get_model_count returns 5") {
    REQUIRE(get_model_count() == 5);
}

TEST_CASE("get_model_names returns 5 names with LinearElastic first") {
    std::vector<std::string> names = get_model_names();
    REQUIRE(names.size() == 5);
    REQUIRE(names[0] == "LinearElastic");
}

TEST_CASE("create_model valid name") {
    // Create a valid MCC model.
    Parameters params(5);
    params << 0.92, 0.2, 1.195, 0.08, 0.02;
    State state(1);
    state << 10.0;

    auto model = create_model("MCC", params, state);
    REQUIRE(model != nullptr);
}

TEST_CASE("create_model invalid name throws") {
    Parameters params = Parameters::Zero(1);
    State state = State::Zero(0);

    REQUIRE_THROWS_AS(create_model("NonExistent", params, state), std::runtime_error);
}

TEST_CASE("create_model_for_metadata by index") {
    // Indices 1 through 5 should succeed.
    for (int i = 1; i <= 5; i++) {
        auto model = create_model_for_metadata(i);
        REQUIRE(model != nullptr);
    }
}

TEST_CASE("create_model_for_metadata out of range throws") {
    REQUIRE_THROWS_AS(create_model_for_metadata(0), std::runtime_error);
    REQUIRE_THROWS_AS(create_model_for_metadata(6), std::runtime_error);
}

TEST_CASE("Each model: param names size matches parameters required") {
    std::vector<std::string> names = get_model_names();
    for (const auto& name : names) {
        auto model = create_model_for_metadata(name);
        REQUIRE(model->get_param_names().size() == static_cast<size_t>(model->get_parameters_required()));
    }
}

TEST_CASE("Each model: state var names size matches state required") {
    std::vector<std::string> names = get_model_names();
    for (const auto& name : names) {
        auto model = create_model_for_metadata(name);
        REQUIRE(model->get_state_var_names().size() == static_cast<size_t>(model->get_state_required()));
    }
}
