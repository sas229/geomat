#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "abstract.hpp"
#include "utilities.hpp"

namespace py = pybind11;

class PyModel : public Model {
    public:
    using Model::Model;

    State get_state_variables() override {
        PYBIND11_OVERRIDE_PURE(
            State,
            Model,
            get_state_variables
        );
    }

    void set_state_variables(
        State new_state
    ) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            Model,
            set_state_variables,
            new_state
        );
    }
};

class PyElastic : public Elastic {
    public:
    using Elastic::Elastic;

    State get_state_variables() override {
        PYBIND11_OVERRIDE_PURE(
            State,
            Elastic,
            get_state_variables
        );
    }

    void set_state_variables(
        State new_state
    ) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            Elastic,
            set_state_variables,
            new_state
        );
    }

    Constitutive compute_D_e(
        Cauchy sigma_prime, 
        Cauchy Delta_epsilon=Cauchy::Zero()
    ) override {
        PYBIND11_OVERRIDE_PURE(
            Constitutive,
            Elastic,
            compute_D_e,
            sigma_prime,
            Delta_epsilon
        );
    }
};

class PyElastoplastic : public Elastoplastic {
    public:
    using Elastoplastic::Elastoplastic;

    State get_state_variables() override {
        PYBIND11_OVERRIDE_PURE(
            State,
            Elastoplastic,
            get_state_variables
        );
    }

    void set_state_variables(
        State new_state
    ) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            Elastoplastic,
            set_state_variables,
            new_state
        );
    }
    
    Constitutive compute_D_e(
        Cauchy sigma_prime, 
        Cauchy Delta_epsilon=Cauchy::Zero()
    ) override {
        PYBIND11_OVERRIDE_PURE(
            Constitutive,
            Elastoplastic,
            compute_D_e,
            sigma_prime,
            Delta_epsilon
        );
    }

    double compute_f(
        Cauchy sigma_prime, 
        State state
    ) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            Elastoplastic,
            compute_f,
            sigma_prime,
            state
        );
    }

    Derivatives compute_derivatives(
        Cauchy sigma_prime, 
        State state
    ) override {
        PYBIND11_OVERRIDE_PURE(
            Derivatives,
            Elastoplastic,
            compute_derivatives,
            sigma_prime,
            state
        );
    }

    State compute_elastic_state_variable_increment(
        Cauchy sigma_prime, 
        State state, 
        Voigt Delta_epsilon_tilde_e
    ) override {
        PYBIND11_OVERRIDE(
            State,
            Elastoplastic,
            compute_elastic_state_variable_increment,
            sigma_prime,
            state,
            Delta_epsilon_tilde_e
        );
    }
};

PYBIND11_MODULE(abstract, m) {
    // Docs.
    m.doc() = "Abstract base classes.";

    // Abstract base classes.
    py::class_<Model, PyModel>(m, "Model")
        .def(py::init<>())
        .def("initialise_log", &Model::initialise_log)
        .def("check_inputs", &Model::check_inputs)
        .def("set_model_name", &Model::set_model_name)
        .def("set_model_type", &Model::set_model_type)
        .def("get_model_name", &Model::get_model_name)
        .def("get_model_type", &Model::get_model_type)
        .def("set_sigma_prime_tilde", &Model::set_sigma_prime_tilde)
        .def("get_sigma_prime_tilde", &Model::get_sigma_prime_tilde)
        .def("set_Delta_epsilon_tilde", &Model::set_Delta_epsilon_tilde)
        .def("get_p_prime", &Model::get_p_prime)
        .def("compute_p_prime", &Model::compute_p_prime)
        .def("compute_q", &Model::compute_q)
        .def("compute_s", &Model::compute_s)
        .def("compute_dq_dsigma_prime", &Model::compute_dq_dsigma_prime)
        .def("get_q", &Model::get_q)
        .def_property("name", &Model::get_model_name, &Model::set_model_name)
        .def_property("type", &Model::get_model_type, &Model::set_model_type)
        .def_property_readonly("p_prime", &Model::get_p_prime)
        .def_property_readonly("q", &Model::get_q)
        .def_property_readonly("sigma_prime_tilde", &Model::get_sigma_prime_tilde)
        .def_property_readonly("I_1", &Model::get_I_1)
        .def_property_readonly("I_2", &Model::get_I_2)
        .def_property_readonly("I_3", &Model::get_I_3)
        .def_property_readonly("J_1", &Model::get_J_1)
        .def_property_readonly("J_2", &Model::get_J_2)
        .def_property_readonly("J_3", &Model::get_J_3)
        .def("solve", &Model::solve);

    py::class_<Elastic, PyElastic, Model>(m, "Elastic", py::multiple_inheritance())
        .def(py::init<>())
        .def("compute_Delta_epsilon_vol", &Elastic::compute_Delta_epsilon_vol)
        .def("compute_K_Butterfield", &Elastic::compute_K_Butterfield)
        .def("compute_G_given_K_nu", &Elastic::compute_G_given_K_nu)
        .def("compute_isotropic_linear_elastic_matrix", &Elastic::compute_isotropic_linear_elastic_matrix)
        .def("get_state_variables", &Elastic::get_state_variables)
        .def("set_state_variables", &Elastic::set_state_variables);

    py::class_<Elastoplastic, PyElastoplastic, Elastic>(m, "Elastoplastic", py::multiple_inheritance())
        .def(py::init<>())
        .def_property("settings", &Elastoplastic::get_settings, &Elastoplastic::set_settings)
        .def("get_state_variables", &Elastoplastic::get_state_variables)
        .def("set_state_variables", &Elastoplastic::set_state_variables);
}