#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <iostream>
#include <fast-forward.hpp>
#include "graph-model-spec.hpp"

namespace py = pybind11;
using namespace clonal_expansion;

//PYBIND11_MAKE_OPAQUE(std::vector<real_t>); // TODO python bindings for
//std::vector<real_t> and a conversion function for the initialiser

PYBIND11_MODULE(pybinding, m) {
    m.doc() = "Pybindings for ff";

    // Bind the Model class
    py::class_<Model>(m, "Model")
        .def(py::init<size_t>(), py::arg("n_vertices"))
        .def_readwrite("m_stages", &Model::m_stages)
        .def_readwrite("m_birth", &Model::m_birth)
        .def_readwrite("m_death", &Model::m_death)
        .def_readwrite("m_immig_rates", &Model::m_immig_rates)
        .def_readwrite("m_initial_pops", &Model::m_initial_pops)
        .def_property("m_migr",
            [](const Model &m) { return m.m_migr; },
            [](Model &m, const std::vector<std::map<int, real_t>> &migr) { m.m_migr = migr; }
        );

    //test binding
    m.def("test_reference", fast_forward::test_reference, "Test", py::arg("qcoords"));

    // Bind rhs_flow
    m.def("rhs_flow", fast_forward::rhs_flow, "Compute rates of change", py::arg("qcoords"), py::arg("parameters"));

    // Bind heun_q_step
    m.def("heun_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::heun_q_step(qcoords, time, dt, parameters);
        for (size_t i = 0; i < qcoords.size(); ++i) {
            std::cout << "qcoords[" << i << "] after Heun step: " << qcoords[i] << std::endl;
        }
    }, "Heun's method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind implicit_q_step
    m.def("implicit_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::implicit_q_step(qcoords, time, dt, parameters);
        for (size_t i = 0; i < qcoords.size(); ++i) {
            std::cout << "qcoords[" << i << "] after Implicit Euler step: " << qcoords[i] << std::endl;
        }
    }, "Implicit Euler method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));
}

