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
    py::class_<C2MC, Elastoplastic, std::shared_ptr<C2MC>>(m, "C2MC", py::multiple_inheritance())
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Extended Mohr Coulomb (EMC).
    py::class_<EMC, Elastoplastic, std::shared_ptr<EMC>>(m, "EMC", py::multiple_inheritance())
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Linear isotropic elasticity (LinearElastic).
    py::class_<LinearElastic, Elastic, std::shared_ptr<LinearElastic>>(m, "LinearElastic", py::multiple_inheritance())
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Modified Cam Clay (MCC).
    py::class_<MCC, Elastoplastic, std::shared_ptr<MCC>>(m, "MCC", py::multiple_inheritance())
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.

    // Soft Modified Cam Clay (SMCC).
    py::class_<SMCC, Elastoplastic, std::shared_ptr<SMCC>>(m, "SMCC", py::multiple_inheritance())
        .def(py::init<Parameters, State>(), py::kw_only(), py::arg("parameters"), py::arg("state")) // Constructor.
        .def(py::init<Parameters, State, std::string>(), py::kw_only(), py::arg("parameters"), py::arg("state"), py::arg("log_severity")); // Overloaded constructor.
}