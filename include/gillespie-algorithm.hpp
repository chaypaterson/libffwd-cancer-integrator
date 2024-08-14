#include <cstddef>
#include <vector>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <graph-model-spec.hpp>
#ifndef GILLESPIE_DEF
#define GILLESPIE_DEF

/* gillespie-algorithm.hpp
 * This is a library that defines an iterator for the Gillespie algorithm when
 * used to simulate birth-death-migration processes on graphs.
 */

namespace clonal_expansion {

namespace gillespie_ssa {

class gillespie_instance {
  public:
    // member variables:
    double m_time{0}; // the simulation time variable
    size_t m_vertices; // the number of vertices on the graph
    Model m_parameters; // the rate parameters and connectivity of the graph
    // the Model class is defined in the header graph-model-spec.hpp
    std::vector<int> m_pops; // the populations on each vertex of the graph
    double m_gamma{0};

    // constructor:
    inline gillespie_instance(const Model &parameters) :
        m_vertices{parameters.m_stages},
        m_parameters{parameters} {
        for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
            // set initial pops:
            m_pops.push_back((int)(parameters.m_initial_pops[vertex]));
        }
        set_gamma(); // FIXME this line causes a segfault in pybind11
    }

    // methods:
    void gillespie_step(gsl_rng *rng); // perform one step of the Gillespie
    // algorithm, updating the time and populations.

  private:
    // These are private methods that are only called from within
    // gillespie_step:
    double get_gamma(void); // return m_gamma
    void set_gamma(void); // compute m_gamma
    std::vector<int> x_to_event(double x); // convert a random variable x to an
    // event (change in populations)
};

double first_passage_time(gsl_rng *rng, const Model &params,
                          const int final_vertex);

std::pair<double,int> first_passage_time(gsl_rng *rng, const Model &params,
        const std::vector<int> final_vertices);

void times_to_final_vertex(const Model &model, int seed,
                           int runs_per_thr, int final_vertex,
                           std::vector<double> &results);

void times_to_final_vertices(const Model &model, int seed,
                             int runs_per_thr, std::vector<int> final_vertices,
                             std::vector<std::pair<double,int>> &results);

void print_results(std::vector<double> &all_times);

void print_kaplan_meier(double time_max, std::vector<double> &all_times,
                        size_t ref_pop);
void print_kaplan_meier(double time_max, std::vector<double> &all_times);

real_t surv_kaplan_meier(double age, std::vector<double> &all_times,
                         size_t ref_pop);

void print_naive_estimator(double time_max, std::vector<double> &all_times);

}

}

#endif //GILLESPIE_DEF
