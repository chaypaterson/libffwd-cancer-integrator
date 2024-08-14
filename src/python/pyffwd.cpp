#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <fast-forward.hpp>
#include <graph-model-spec.hpp>
#include <gillespie-algorithm.hpp>
#include <gsl/gsl_rng.h>

namespace py = pybind11;
using namespace clonal_expansion;

PYBIND11_MAKE_OPAQUE(std::vector<real_t>);

// GSL_RNG factory function:
gsl_rng* seed_gsl_rng(int seed) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    return r;
}

PYBIND11_MODULE(pyffwd, m) {
    m.doc() = "Python bindings for the cancer-integrator (fastforward/ffwd) library";

    // Convert Python list to std::vector<real_t>
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

    // Bind the rhs_flow function
    m.def("rhs_flow", fast_forward::rhs_flow, 
          "Compute rates of change", py::arg("qcoords"), py::arg("parameters"));

    // Bind the heun_q_step function
    m.def("heun_q_step", fast_forward::heun_q_step, 
          "Heun's method for q-coordinates", 
          py::arg("qcoords"), py::arg("time"), py::arg("dt"),
          py::arg("parameters"));

    // Bind the implicit_q_step function
    m.def("implicit_q_step", fast_forward::implicit_q_step, 
          "Implicit Euler method for q-coordinates", 
          py::arg("qcoords"), py::arg("time"), py::arg("dt"),
          py::arg("parameters"));

    // Bind the rungekutta_q_step function
    m.def("rungekutta_q_step", fast_forward::rungekutta_q_step, 
          "Runge-Kutta method for q-coordinates", 
          py::arg("qcoords"), py::arg("time"), py::arg("dt"),
          py::arg("parameters"));

    // Bind the generating_function function
    m.def("generating_function", fast_forward::generating_function, 
          "Compute the generating function", 
          py::arg("qcoords"), py::arg("initial_pops"));

    // Bind the gillespie_instance class
    using gillespie_ssa::gillespie_instance;
    py::class_<gillespie_instance>(m, "GillespieInstance")
        .def(py::init<const Model &>(), py::arg("model"))
        .def("gillespie_step", &gillespie_instance::gillespie_step)
        .def_readwrite("m_time", &gillespie_instance::m_time)
        .def_readwrite("m_vertices", &gillespie_instance::m_vertices)
        .def_readwrite("m_pops", &gillespie_instance::m_pops)
        .def_readwrite("m_parameters", &gillespie_instance::m_parameters);

    // Bind gsl_rng as a class
    py::class_<gsl_rng>(m, "GSL_RNG")
        .def(py::init<>())
        .def_readwrite("type", &gsl_rng::type)
        .def_readwrite("state", &gsl_rng::state);

    // Bind gsl_rng helper functions:
    m.def("seed_gsl_rng", seed_gsl_rng, 
          "Create new GSL RNG object",
          py::arg("seed")); // TODO: make into custom constructor for GSL_RNG?
    //m.def("free_gsl_rng", gsl_rng_free,
    //      "Free GSL RNG object");

    // Bind standalone functions
    m.def("first_passage_time_multiple", 
          static_cast<std::pair<double,int> (*)(gsl_rng *, const Model&,
                                                const std::vector<int>)>(&gillespie_ssa::first_passage_time));

    m.def("times_to_final_vertex", gillespie_ssa::times_to_final_vertex, 
          "Compute times to reach the final vertex across multiple runs", 
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"),
          py::arg("final_vertex"), py::arg("results"));

    m.def("times_to_final_vertices", &gillespie_ssa::times_to_final_vertices,
          "Compute times to reach any of the final vertices across multiple runs",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"),
          py::arg("final_vertices"), py::arg("results"));

    // Bind print_results function
    m.def("print_results", [](const std::vector<double> &all_times) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        clonal_expansion::gillespie_ssa::print_results(non_const_times);
    }, "Print results");

    // Bind print_kaplan_meier
    m.def("print_kaplan_meier", [](double time_max, const std::vector<double> &all_times, size_t ref_pop) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        clonal_expansion::gillespie_ssa::print_kaplan_meier(time_max, non_const_times, ref_pop);
    }, "Print Kaplan-Meier with ref_pop", py::arg("time_max"), py::arg("all_times"), py::arg("ref_pop"));

    // Bind surv_kaplan_meier
    m.def("surv_kaplan_meier", [](double age, const std::vector<double> &all_times, size_t ref_pop) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        return clonal_expansion::gillespie_ssa::surv_kaplan_meier(age, non_const_times, ref_pop);
    }, "Survival Kaplan-Meier", py::arg("age"), py::arg("all_times"), py::arg("ref_pop"));
    
    // Bind print_naive_estimator
    m.def("print_naive_estimator", [](double time_max, const std::vector<double> &all_times) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        clonal_expansion::gillespie_ssa::print_naive_estimator(time_max, non_const_times);
    }, "Print naive estimator", py::arg("time_max"), py::arg("all_times"));

}
