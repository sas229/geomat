#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Eigen/Eigen>

#include "utilities.hpp"

namespace py = pybind11;

PYBIND11_MODULE(utilities, m) {
    // Docs.
    m.doc() = "Model implementation classes.";

    // Settings.
    py::class_<Settings>(m, "Settings")
        .def(py::init<>())
        .def_property("FTOL", &Settings::get_FTOL, &Settings::set_FTOL)
        .def_property("LTOL", &Settings::get_LTOL, &Settings::set_LTOL)
        .def_property("STOL", &Settings::get_STOL, &Settings::set_STOL)
        .def_property("EPS", &Settings::get_EPS, &Settings::set_EPS)
        .def_property("DT_MIN", &Settings::get_DT_MIN, &Settings::set_DT_MIN)
        .def_property("MAXITS_YSI", &Settings::get_MAXITS_YSI, &Settings::set_MAXITS_YSI)
        .def_property("MAXITS_YSC", &Settings::get_MAXITS_YSC, &Settings::set_MAXITS_YSC)
        .def_property("NSUB", &Settings::get_NSUB, &Settings::set_NSUB);

    // Derivatives object.
    py::class_<Derivatives>(m, "Derivatives")
        .def(py::init<>())
        .def_property("df_dsigma_prime", &Derivatives::get_df_dsigma_prime, &Derivatives::set_df_dsigma_prime)
        .def_property("dg_dsigma_prime", &Derivatives::get_dg_dsigma_prime, &Derivatives::set_dg_dsigma_prime)
        .def_property("H_s", &Derivatives::get_H_s, &Derivatives::set_H_s)
        .def_property("B_s", &Derivatives::get_B_s, &Derivatives::set_B_s);
}