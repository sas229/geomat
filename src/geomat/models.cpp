#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "Model.hpp"
#include "Elastic.hpp"
#include "Elastoplastic.hpp"

#include "LinearElastic.hpp"
#include "MCC.hpp"
#include "SMCC.hpp"

namespace py = pybind11;

PYBIND11_MODULE(models, m) {
    // // Docs.
    // m.doc() = "Constitutive model library.";

    // // Standalone module definition.
    // auto package = pybind11::module::import("package");
    // auto module =  package.attr("module");
    // m.add_object("module", module);

    // Abstract base classes.
    py::class_<Model>(m, "Model");
    py::class_<Elastic, Model>(m, "Elastic");
    py::class_<Elastoplastic, Elastic>(m, "Elastoplastic");

    // Models.

    // Linear isotropic elasticity (LinearElastic).
    py::class_<LinearElastic, Elastic>(m, "LinearElastic")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")) // Overloaded constructor.
        .def("set_sigma_prime_tilde", &LinearElastic::set_sigma_prime_tilde)
        .def("get_sigma_prime_tilde", &LinearElastic::get_sigma_prime_tilde)
        .def("set_Delta_epsilon_tilde", &LinearElastic::set_Delta_epsilon_tilde)
        .def("get_state", &LinearElastic::get_state_variables)
        .def("get_p_prime", &LinearElastic::get_p_prime)
        .def("get_q", &LinearElastic::get_q)
        .def_property_readonly("name", &LinearElastic::get_model_name)
        .def_property_readonly("type", &LinearElastic::get_model_type)
        .def_property_readonly("p_prime", &LinearElastic::get_p_prime)
        .def_property_readonly("q", &LinearElastic::get_q)
        .def_property_readonly("sigma_prime_tilde", &LinearElastic::get_sigma_prime_tilde)
        .def_property_readonly("state", &LinearElastic::get_state_variables)
        .def_property_readonly("I_1", &LinearElastic::get_I_1)
        .def_property_readonly("I_2", &LinearElastic::get_I_2)
        .def_property_readonly("I_3", &LinearElastic::get_I_3)
        .def_property_readonly("J_1", &LinearElastic::get_J_1)
        .def_property_readonly("J_2", &LinearElastic::get_J_2)
        .def_property_readonly("J_3", &LinearElastic::get_J_3)
        .def("solve", &LinearElastic::solve);

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic>(m, "MCC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")) // Overloaded constructor.
        .def("set_sigma_prime_tilde", &MCC::set_sigma_prime_tilde)
        .def("get_sigma_prime_tilde", &MCC::get_sigma_prime_tilde)
        .def("set_Delta_epsilon_tilde", &MCC::set_Delta_epsilon_tilde)
        .def("get_state", &MCC::get_state_variables)
        .def("get_p_prime", &MCC::get_p_prime)
        .def("get_q", &MCC::get_q)
        .def_property_readonly("name", &MCC::get_model_name)
        .def_property_readonly("type", &MCC::get_model_type)
        .def_property_readonly("p_prime", &MCC::get_p_prime)
        .def_property_readonly("q", &MCC::get_q)
        .def_property_readonly("sigma_prime_tilde", &MCC::get_sigma_prime_tilde)
        .def_property_readonly("state", &MCC::get_state_variables)
        .def_property_readonly("I_1", &MCC::get_I_1)
        .def_property_readonly("I_2", &MCC::get_I_2)
        .def_property_readonly("I_3", &MCC::get_I_3)
        .def_property_readonly("J_1", &MCC::get_J_1)
        .def_property_readonly("J_2", &MCC::get_J_2)
        .def_property_readonly("J_3", &MCC::get_J_3)
        .def("solve", &MCC::solve);

    // Soft Modified Cam Clay (SMCC).
    py::class_<SMCC, Elastoplastic>(m, "SMCC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")) // Overloaded constructor.
        .def("set_sigma_prime_tilde", &SMCC::set_sigma_prime_tilde)
        .def("get_sigma_prime_tilde", &SMCC::get_sigma_prime_tilde)
        .def("set_Delta_epsilon_tilde", &SMCC::set_Delta_epsilon_tilde)
        .def("get_state", &SMCC::get_state_variables)
        .def("get_p_prime", &SMCC::get_p_prime)
        .def("get_q", &SMCC::get_q)
        .def_property_readonly("name", &SMCC::get_model_name)
        .def_property_readonly("type", &SMCC::get_model_type)
        .def_property_readonly("p_prime", &SMCC::get_p_prime)
        .def_property_readonly("q", &SMCC::get_q)
        .def_property_readonly("sigma_prime_tilde", &SMCC::get_sigma_prime_tilde)
        .def_property_readonly("state", &SMCC::get_state_variables)
        .def_property_readonly("I_1", &SMCC::get_I_1)
        .def_property_readonly("I_2", &SMCC::get_I_2)
        .def_property_readonly("I_3", &SMCC::get_I_3)
        .def_property_readonly("J_1", &SMCC::get_J_1)
        .def_property_readonly("J_2", &SMCC::get_J_2)
        .def_property_readonly("J_3", &SMCC::get_J_3)
        .def("solve", &SMCC::solve);
}