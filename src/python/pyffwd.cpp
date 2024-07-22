#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <iostream>
#include <fast-forward.hpp>
#include <graph-model-spec.hpp>

namespace py = pybind11;
using namespace clonal_expansion;
using clonal_expansion::real_t;

PYBIND11_MAKE_OPAQUE(std::vector<real_t>);


PYBIND11_MODULE(pyffwd, m) {
    m.doc() = "Pybindings for ff";

    // convert Python list to std::vector<real_t>
    m.def("list_to_vector", [](py::list py_list) {
        std::vector<real_t> cpp_vector;
        for (auto item : py_list) {
            cpp_vector.push_back(item.cast<real_t>());
        }
        return cpp_vector;
    }, "Convert list to std::vector<real_t>", py::arg("py_list"));
    
    py::bind_vector<std::vector<real_t>>(m, "RealVector");

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

    // Bind rhs_flow
    m.def("rhs_flow", [](const std::vector<real_t> &qcoords, Model &parameters) {
    
        auto result = fast_forward::rhs_flow(qcoords, parameters);
        return result;
    }, "Compute rates of change", py::arg("qcoords"), py::arg("parameters"));

    // Bind heun_q_step
    m.def("heun_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::heun_q_step(qcoords, time, dt, parameters);
    }, "Heun's method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind implicit_q_step
    m.def("implicit_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::implicit_q_step(qcoords, time, dt, parameters);
    }, "Implicit Euler method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind rungekutta_q_step
    m.def("rungekutta_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::rungekutta_q_step(qcoords, time, dt, parameters);
    }, "Runge-Kutta method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind generating_function
    m.def("generating_function", [](const std::vector<real_t> &qcoords, const std::vector<real_t> &initial_pops) {
        return fast_forward::generating_function(qcoords, initial_pops);
    }, "Compute the generating function", py::arg("qcoords"), py::arg("initial_pops"));

}

