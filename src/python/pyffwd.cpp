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
using clonal_expansion::real_t;

PYBIND11_MAKE_OPAQUE(std::vector<real_t>);

// Define RNG wrapper class
class RNGWrapper {
public:
    RNGWrapper(unsigned long seed) {
        rng = gsl_rng_alloc(gsl_rng_default); //allocate a random number generator (RNG)
        gsl_rng_set(rng, seed); //The seed is used to initialize the RNG
    }

    ~RNGWrapper() {
        gsl_rng_free(rng); //Releases the memory allocated for the RNG using GSL
    }

    // calculate the first passage time to a single vertex
    //return gillespie_ssa::first_passage_time(rng, model, final_vertex)
    double first_passage_time_single(const Model& model, int final_vertex) {
        return gillespie_ssa::first_passage_time(rng, model, final_vertex);
    }

    //calculate the first passage times for multiple final vertices
    std::vector<std::pair<double, int>> first_passage_time_multiple(const Model& model, const std::vector<int>& final_vertices) {
        std::vector<std::pair<double, int>> results;
        for (int final_vertex : final_vertices) {
            std::pair<double, int> result = gillespie_ssa::first_passage_time(rng, model, final_vertices);
            if (result.first >= 0) { // Only add valid results
                results.push_back(result);
            }
        }
        return results;
    }

private:
    gsl_rng* rng; //Declares a pointer to a gsl_rng object. This pointer holds the reference to the RNG allocated by GSL
};

PYBIND11_MODULE(pyffwd, m) {
    m.doc() = "Pybindings for ff";

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
    m.def("rhs_flow", [](const std::vector<real_t> &qcoords, Model &parameters) {
        auto result = fast_forward::rhs_flow(qcoords, parameters);
        return result;
    }, "Compute rates of change", py::arg("qcoords"), py::arg("parameters"));

    // Bind the heun_q_step function
    m.def("heun_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::heun_q_step(qcoords, time, dt, parameters);
    }, "Heun's method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind the implicit_q_step function
    m.def("implicit_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::implicit_q_step(qcoords, time, dt, parameters);
    }, "Implicit Euler method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind the rungekutta_q_step function
    m.def("rungekutta_q_step", [](std::vector<real_t> &qcoords, const real_t &time, real_t &dt, Model &parameters) {
        fast_forward::rungekutta_q_step(qcoords, time, dt, parameters);
    }, "Runge-Kutta method for q-coordinates", py::arg("qcoords"), py::arg("time"), py::arg("dt"), py::arg("parameters"));

    // Bind the generating_function function
    m.def("generating_function", [](const std::vector<real_t> &qcoords, const std::vector<real_t> &initial_pops) {
        return fast_forward::generating_function(qcoords, initial_pops);
    }, "Compute the generating function", py::arg("qcoords"), py::arg("initial_pops"));

    // Bind the gillespie_instance class
    py::class_<gillespie_ssa::gillespie_instance>(m, "GillespieInstance")
        .def(py::init<const Model &>(), py::arg("model"))
        .def("gillespie_step", &gillespie_ssa::gillespie_instance::gillespie_step, py::arg("rng"))
        .def_readwrite("m_pops", &gillespie_ssa::gillespie_instance::m_pops)
        .def_readwrite("m_time", &gillespie_ssa::gillespie_instance::m_time);

    // Bind the RNGWrapper class
    py::class_<RNGWrapper>(m, "RNGWrapper")
        .def(py::init<unsigned long>(), py::arg("seed"))
        .def("first_passage_time_single", &RNGWrapper::first_passage_time_single, py::arg("model"), py::arg("final_vertex"))
        .def("first_passage_time_multiple", &RNGWrapper::first_passage_time_multiple, py::arg("model"), py::arg("final_vertices"));

    // Bind standalone functions
    m.def("times_to_final_vertex", &gillespie_ssa::times_to_final_vertex,
          "Compute times to reach the final vertex across multiple runs",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"), py::arg("final_vertex"), py::arg("results"));

    m.def("times_to_final_vertices", &gillespie_ssa::times_to_final_vertices,
          "Compute times to reach any of the final vertices across multiple runs",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"), py::arg("final_vertices"), py::arg("results"));

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
