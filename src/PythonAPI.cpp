#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "Model.hpp"
#include "Elastic.hpp"
#include "Elastoplastic.hpp"

#include "LinearElastic.hpp"
#include "MCC.hpp"

namespace py = pybind11;

PYBIND11_MODULE(models, m) {
    // Abstract base classes.
    py::class_<Model>(m, "Model");
    py::class_<Elastic, Model>(m, "Elastic");
    py::class_<Elastoplastic, Elastic>(m, "Elastoplastic");

    // Models.

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic>(m, "MCC")
        .def(py::init<Parameters, State>()) // Constructor.
        .def("set_sigma_prime_tilde", &MCC::set_sigma_prime_tilde)
        .def("get_sigma_prime_tilde", &MCC::get_sigma_prime_tilde)
        .def("get_p_prime", &MCC::get_p_prime)
        .def("get_q", &MCC::get_q)
        .def_property_readonly("p_prime", &MCC::get_p_prime)
        .def_property_readonly("q", &MCC::get_q)
        .def_property_readonly("sigma_prime_tilde", &MCC::get_sigma_prime_tilde)
        .def_property_readonly("I_1", &MCC::get_I_1)
        .def_property_readonly("I_2", &MCC::get_I_2)
        .def_property_readonly("I_3", &MCC::get_I_3)
        .def_property_readonly("J_1", &MCC::get_J_1)
        .def_property_readonly("J_2", &MCC::get_J_2)
        .def_property_readonly("J_3", &MCC::get_J_3)
        .def("set_Delta_epsilon_tilde", &MCC::set_Delta_epsilon_tilde)
        .def("solve", &MCC::solve);
}