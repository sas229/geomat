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

    // Isotropic linear elasticity.
    py::class_<LinearElastic, Elastic>(m, "LinearElastic");
    //     m.def(py::init<Parameters, State>()); // Constructor.

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic>(m, "MCC")
        .def(py::init<Parameters, State>()) // Constructor.
        .def("set_sigma_prime", &MCC::set_sigma_prime)
        .def("get_sigma_prime", &MCC::get_sigma_prime)
        .def("set_strain_increment", &MCC::set_strain_increment)
        .def("solve", &MCC::solve);
}