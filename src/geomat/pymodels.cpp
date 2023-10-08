#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "abstract.hpp"
#include "models.hpp"

namespace py = pybind11;

PYBIND11_MODULE(models, m) {
    // Docs.
    m.doc() = "Model implementation classes.";

    // C2 Continuous Mohr Coulomb (C2MC).
    py::class_<C2MC, Elastoplastic>(m, "C2MC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Extended Mohr Coulomb (EMC).
    py::class_<EMC, Elastoplastic>(m, "EMC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Linear isotropic elasticity (LinearElastic).
    py::class_<LinearElastic, Elastic>(m, "LinearElastic")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic>(m, "MCC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Soft Modified Cam Clay (SMCC).
    py::class_<SMCC, Elastoplastic>(m, "SMCC")
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.
}