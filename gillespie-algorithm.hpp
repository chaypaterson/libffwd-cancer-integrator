#include <cstddef>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "graph-model-spec.hpp"
#ifndef GILLESPIE_DEF
#define GILLESPIE_DEF

/* gillespie-algorithm.hpp
 * This is a library that defines an iterator for the Gillespie algorithm when
 * used to simulate birth-death-migration processes on graphs.
 */
class gillespie_instance {
public:
    // member variables:
    double m_time{0}; // the simulation time variable
    std::vector<int> m_pops; // the populations on each vertex of the graph
    Model m_parameters; // the rate parameters and connectivity of the graph
    // the Model class is defined in the header graph-model-spec.hpp
    size_t m_vertices; // the number of vertices on the graph

    // constructor:
    inline gillespie_instance(const Model &parameters) :
        m_vertices{parameters.m_stages},
        m_parameters{parameters} {
        for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
            // set initial pops:
            m_pops[vertex] = (int)parameters.m_initial_pops[vertex];
        }
    }
    
    // methods:
    void gillespie_step(gsl_rng *rng); // perform one step of the Gillespie
    // algorithm, updating the time and populations.

private:
    // These are private methods that are only called from within
    // gillespie_step:
    double get_gamma(void); // compute Gamma
    std::vector<int> x_to_event(double x); // convert a random variable x to an
    // event (change in populations)
};

#endif
