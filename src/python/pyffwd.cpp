#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <fast-forward.hpp>
#include <graph-model-spec.hpp>
#include <gillespie-algorithm.hpp>
#include <gsl/gsl_rng.h>
#include <pybind11/pytypes.h>

namespace py = pybind11;
using namespace clonal_expansion;

PYBIND11_MAKE_OPAQUE(std::vector<real_t>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);

// GSL_RNG factory function:
gsl_rng* seed_gsl_rng(int seed) {
    const gsl_rng_type *T = gsl_rng_mt19937;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);


    return r;
}

// Convert std::vector<int> to Python list
py::list vector_to_pylist(const std::vector<int>& vec) {
    py::list pylist;
    for (const int& item : vec) {
        pylist.append(item);
    }
    return pylist;
}

// Convert Python list to std::vector<int>
std::vector<int> pylist_to_vector(const py::list& pylist) {
    std::vector<int> vec;
    for (const auto& item : pylist) {
        vec.push_back(item.cast<int>());
    }
    return vec;
}

// Convert Python list to std::vector<real_t>
std::vector<real_t> vec_real_t(const py::list& pylist) {
    std::vector<real_t> vec;
    for (const auto& item : pylist) {
        vec.push_back(item.cast<real_t>());
    }
    return vec;
}

PYBIND11_MODULE(pyffwd, m) {
    m.doc() = "Python bindings for the cancer-integrator (fastforward/ffwd) library";

    // Convert Python list to std::vector<real_t>
    m.def("list_to_vector", [](py::list py_list) {
        return vec_real_t(py_list);
    }, "Convert list to std::vector<real_t>", py::arg("py_list"));

    py::bind_vector<std::vector<real_t>>(m, "RealVector");

    // Convert Python list to std::vector<int>
    m.def("list_vector_int", [](py::list py_list) {
        std::vector<int> cpp_vector;
        for (auto item : py_list) {
            cpp_vector.push_back(item.cast<int>());
        }
        return cpp_vector;
    }, "Convert list to std::vector<int>", py::arg("py_list"));

    py::bind_vector<std::vector<int>>(m, "IntVector");

    // Bind the Model class
    py::class_<Model>(m, "Model")
         // Original constructor
        .def(py::init<size_t>(), py::arg("n_vertices"))
        // Factory function for pass birth, death rates, initial_pops
        .def_static("create_model",
            [](size_t n_vertices,
               py::object birth_rates = py::none(),
               py::object death_rates = py::none(),
               py::object initial_pops = py::none()) {

                // Create a default Model
                Model model(n_vertices);

                if (!birth_rates.is_none()) {
                    auto birth_list = vec_real_t(birth_rates.cast<py::list>());
                    model.m_birth = birth_list;
                }

                if (!death_rates.is_none()) {
                    auto death_list = vec_real_t(death_rates.cast<py::list>());
                    model.m_death = death_list;
                }

                if (!initial_pops.is_none()) {
                    auto pops_list = vec_real_t(initial_pops.cast<py::list>());
                    model.m_initial_pops = pops_list;
                }

                return model;
            },
            py::arg("n_vertices"),
            py::arg("birth_rates") = py::none(),
            py::arg("death_rates") = py::none(),
            py::arg("initial_pops") = py::none())

        // Setter for birth, death, and initial populations
        .def("set_birth", [](Model &model, py::list birth_rates) {
            model.m_birth = vec_real_t(birth_rates);
        })
        .def("set_death", [](Model &model, py::list death_rates) {
            model.m_death = vec_real_t(death_rates);
        })
        .def("set_initial_pops", [](Model &model, py::list initial_pops) {
            model.m_initial_pops = vec_real_t(initial_pops);
        })
        
        // Accessors for internal variables
        .def_readwrite("m_stages", &Model::m_stages)
        .def_readwrite("m_birth", &Model::m_birth)
        .def_readwrite("m_death", &Model::m_death)
        .def_readwrite("m_immig_rates", &Model::m_immig_rates)
        .def_readwrite("m_initial_pops", &Model::m_initial_pops)
        .def_property("m_migr",
            [](const Model &m) { return m.m_migr; },
            [](Model &m, const std::vector<std::map<int, real_t>> &migr) {
                m.m_migr = migr;
            }
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
        .def(py::init([](const Model &params){
                gillespie_instance instance(params);
                return instance;
            }),
             py::arg("model"))
        .def("gillespie_step", &gillespie_instance::gillespie_step)
        .def_readwrite("m_time", &gillespie_instance::m_time)
        .def_readwrite("m_vertices", &gillespie_instance::m_vertices)
        .def_property("m_pops",
            [](const gillespie_instance &instance) {
                return vector_to_pylist(instance.m_pops);
            },
            [](gillespie_instance &instance, const py::list &pylist) {
                instance.m_pops = pylist_to_vector(pylist);
            }
        )
        .def_readwrite("m_parameters", &gillespie_instance::m_parameters);

    // Bind gsl_rng as a class
    py::class_<gsl_rng>(m, "GSL_RNG")
        .def(py::init([](int seed) {
            return seed_gsl_rng(seed);
        }), "Create new GSL RNG object", py::arg("seed"))
        .def("uniform", [](gsl_rng* rng) {
            return gsl_rng_uniform(rng);
        }, "Generate a random number in the interval [0, 1)")
        .def("uniform_int", [](gsl_rng* rng, int max) {
            return gsl_rng_uniform_int(rng, max);
        }, "Generate a random integer in the interval [0, max)")
        .def_readwrite("type", &gsl_rng::type)
        .def_readwrite("state", &gsl_rng::state);

    // Bind standalone functions

    m.def("first_passage_time",
        static_cast<double (*)(gsl_rng *, const Model&, const int)>(
            &gillespie_ssa::first_passage_time),
        py::arg("rng"), py::arg("model"), py::arg("final_vertex"));


    m.def("first_passage_time_multiple",
        // Cast the corresponding function to a function pointer type:
        static_cast<std::pair<double,int> (*)(
            gsl_rng *, const Model&, const std::vector<int>)>(
                &gillespie_ssa::first_passage_time));

    m.def("times_to_final_vertex", gillespie_ssa::times_to_final_vertex,
          "Compute times to reach the final vertex across multiple runs",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"),
          py::arg("final_vertex"), py::arg("results"));

    m.def("times_to_final_vertices", &gillespie_ssa::times_to_final_vertices,
          "Compute times to reach any of the final vertices across multiple runs",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"),
          py::arg("final_vertices"), py::arg("results"));

    m.def("first_passage_time_poly", &gillespie_ssa::first_passage_time_poly,
          "Compute the first passage time to reach a set of final vertices", 
          py::arg("rng"), py::arg("model"), py::arg("final_vertices"));

    m.def("times_to_final_vertices_poly", &gillespie_ssa::times_to_final_vertices_poly,
          "Compute times to reach any of the final vertices across multiple runs (poly)",
          py::arg("model"), py::arg("seed"), py::arg("runs_per_thr"),
          py::arg("final_vertices"), py::arg("results"));

    // Bind print_results function
    m.def("print_results", [](const std::vector<double> &all_times) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        clonal_expansion::gillespie_ssa::print_results(non_const_times);
    }, "Print results");


    // Binding Kaplan-Meier functions
    m.def("print_kaplan_meier",
        // Cast the corresponding function to a function pointer type:
        static_cast<void (*)(double, std::vector<double>&)>(
            &clonal_expansion::gillespie_ssa::print_kaplan_meier),
        py::arg("time_max"), py::arg("all_times"));

    m.def("print_kaplan_meier", 
          // Cast the corresponding function to a function pointer type:
          static_cast<void(*)(double, std::vector<double>&, size_t)>(
              &clonal_expansion::gillespie_ssa::print_kaplan_meier),
          py::arg("time_max"), py::arg("all_times"), py::arg("ref_pop"));

    m.def("surv_kaplan_meier", &clonal_expansion::gillespie_ssa::surv_kaplan_meier, 
          py::arg("age"), py::arg("all_times"), py::arg("ref_pop"));

    // Bind print_naive_estimator
    m.def("print_naive_estimator", [](double time_max,
                                      const std::vector<double> &all_times) {
        auto& non_const_times = const_cast<std::vector<double>&>(all_times);
        clonal_expansion::gillespie_ssa::print_naive_estimator(time_max,
                                                               non_const_times);
    }, "Print naive estimator", py::arg("time_max"), py::arg("all_times"));

}
