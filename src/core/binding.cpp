#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "fast-forward.hpp"       
#include "graph-model-spec.hpp"       

namespace py = pybind11;
using namespace clonal_expansion;

PYBIND11_MODULE(pybinding, m) {
    m.doc() = "Pybindings for ff";

    // Bind rhs_flow
    m.def("rhs_flow", [](const std::vector<real_t> &qcoords, Model &parameters) {
        return fast_forward::rhs_flow(qcoords, parameters);
    }, "Compute rates of change", py::arg("qcoords"), py::arg("parameters"));

    // Bind the Model class 
    py::class_<Model>(m, "Model")
        .def(py::init<size_t>(), py::arg("n_vertices"))
        .def_readwrite("m_stages", &Model::m_stages)
        .def_readwrite("m_birth", &Model::m_birth)
        .def_readwrite("m_death", &Model::m_death)
        .def_readwrite("m_immig_rates", &Model::m_immig_rates)
        .def_readwrite("m_initial_pops", &Model::m_initial_pops)
        // getter and setter for m_migr
        .def_property("m_migr", 
            [](const Model &m) { return m.m_migr; },
            [](Model &m, const std::vector<std::map<int, real_t>> &migr) { m.m_migr = migr; }
        );

    // Bind heun_q_step
    m.def("heun_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::heun_q_step(qcoords, time, dt, parameters);
    }, "Heun's method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));
}
