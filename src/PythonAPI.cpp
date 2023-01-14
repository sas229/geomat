#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "Model.hpp"
#include "Elastic.hpp"
#include "Elastoplastic.hpp"

#include "LinearElastic.hpp"
#include "MCC.hpp"
// #include "SMCC.hpp"

namespace py = pybind11;

PYBIND11_MODULE(models, m) {
    // Abstract base classes.
    py::class_<Model>(m, "Model");
    py::class_<Elastic, Model>(m, "Elastic");
    py::class_<Elastoplastic, Elastic>(m, "Elastoplastic");

    // Models.

    // Isotropic linear elasticity.
    py::class_<LinearElastic, Elastic>(m, "LinearElastic")
        .def(py::init<Parameters, State>()); // Constructor.

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic>(m, "MCC")
        .def(py::init<Parameters, State>()) // Constructor.
        .def(py::init<Parameters, State, std::string>()) // Overloaded constructor.
        .def("set_sigma_prime_tilde", &MCC::set_sigma_prime_tilde)
        .def("set_Delta_epsilon_tilde", &MCC::set_Delta_epsilon_tilde)
        .def("set_IP_number", &MCC::set_IP_number)
        .def("solve", &MCC::solve)
        .def_property_readonly("sigma_prime_tilde", &MCC::get_sigma_prime_tilde)
        .def_property_readonly("sigma_prime", &MCC::get_sigma_prime)
        .def_property_readonly("name", &MCC::get_name)
        .def_property_readonly("model_type", &MCC::get_model_type)
        .def_property_readonly("IP_number", &MCC::get_IP_number)
        .def_property_readonly("p_prime", &MCC::get_p_prime)
        .def_property_readonly("q", &MCC::get_q)
        .def_property_readonly("jacobian", &MCC::get_jacobian)
        .def_property_readonly("I_1", &MCC::get_I_1)
        .def_property_readonly("I_2", &MCC::get_I_2)
        .def_property_readonly("I_3", &MCC::get_I_3)
        .def_property_readonly("J_1", &MCC::get_J_1)
        .def_property_readonly("J_2", &MCC::get_J_2)
        .def_property_readonly("J_3", &MCC::get_J_3)
        .def_property_readonly("solved", &MCC::get_solved)
        .def_property_readonly("mises_stress", &MCC::get_mises_stress)
        .def_property_readonly("max_shear", &MCC::get_max_shear)
        ;
}