#include "ModelLookup.hpp"
#include "models.hpp"
#include <stdexcept>

struct ModelEntry {
    std::string name;
    int n_params;
    int n_state;
};

static const std::vector<ModelEntry> model_entries = {
    {"LinearElastic", 2, 0},
    {"C2MC", 7, 0},
    {"EMC", 8, 1},
    {"MCC", 5, 1},
    {"SMCC", 8, 2},
};

static std::unique_ptr<Model> construct(const std::string& name, Parameters parameters, State state, std::string log_severity) {
    if (name == "LinearElastic") {
        return std::make_unique<LinearElastic>(parameters, state, log_severity);
    } else if (name == "C2MC") {
        return std::make_unique<C2MC>(parameters, state, log_severity);
    } else if (name == "EMC") {
        return std::make_unique<EMC>(parameters, state, log_severity);
    } else if (name == "MCC") {
        return std::make_unique<MCC>(parameters, state, log_severity);
    } else if (name == "SMCC") {
        return std::make_unique<SMCC>(parameters, state, log_severity);
    }
    throw std::runtime_error("Unknown model name: " + name);
}

std::unique_ptr<Model> create_model(const std::string& name, Parameters parameters, State state, std::string log_severity) {
    return construct(name, parameters, state, log_severity);
}

int get_model_count() {
    return static_cast<int>(model_entries.size());
}

std::vector<std::string> get_model_names() {
    std::vector<std::string> names;
    names.reserve(model_entries.size());
    for (const auto& entry : model_entries) {
        names.push_back(entry.name);
    }
    return names;
}

std::unique_ptr<Model> create_model_for_metadata(const std::string& name) {
    for (const auto& entry : model_entries) {
        if (entry.name == name) {
            Parameters params = Parameters::Zero(entry.n_params);
            State state = State::Zero(entry.n_state);
            return construct(name, params, state, "none");
        }
    }
    throw std::runtime_error("Unknown model name: " + name);
}

std::unique_ptr<Model> create_model_for_metadata(int index) {
    if (index < 1 || index > static_cast<int>(model_entries.size())) {
        throw std::runtime_error("Model index out of range: " + std::to_string(index));
    }
    return create_model_for_metadata(model_entries[index - 1].name);
}
