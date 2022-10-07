#include <cstddef>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "graph-model-spec.hpp"
#include "gillespie-algorithm.hpp"

/* Methods for the gillespie_instance class, which represents a dynamical state
 * of a birth-death process on a graph, as simulated by the Gillespie
 * algorithm.
 */

// a method to perform one step, changing the populations and advancing time
void gillespie_instance::gillespie_step(gsl_rng *rng) {
    double Gamma = get_gamma();
    double x = gsl_ran_flat(rng, 0, Gamma);
    auto delta_pops = x_to_event(x);

    for (int vertex = 0; vertex < m_vertices; ++vertex) {
        m_pops[vertex] += delta_pops[vertex];
    }

    double Delta_t = gsl_ran_exponential(rng, 1.0 / Gamma);

    m_time += Delta_t;
}

// a method to compute Gamma
double gillespie_instance::get_gamma() {
    double Gamma = 0;
    // for each possible process, increment Gamma by the event rate = rate
    // parameter * population on this node
    for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
        Gamma += m_parameters.m_birth[vertex] * m_pops[vertex];
        Gamma += m_parameters.m_death[vertex] * m_pops[vertex];
        for (size_t out_vertex = 0; 
             out_vertex < m_vertices; ++out_vertex) {
            Gamma += m_parameters.m_migr[vertex][out_vertex] * m_pops[vertex];
        }
    }
    return Gamma;
}

// a method to accept a variable x < Gamma and return an event, which
// consists of a change in the populations
std::vector<int> gillespie_instance::x_to_event(double x) {
    double accum_rate = 0;
    std::vector<int> delta_pops(m_vertices);
    for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
        accum_rate += m_parameters.m_birth[vertex] * m_pops[vertex];
        if (accum_rate >= x) {
            delta_pops[vertex] = +1;
            break;
        }
        accum_rate += m_parameters.m_death[vertex] * m_pops[vertex];
        if (accum_rate >= x) {
            delta_pops[vertex] = -1;
            break;
        }
        for (size_t out_vertex = 0; 
             out_vertex < m_vertices; ++out_vertex) {
            accum_rate += m_parameters.m_migr[vertex][out_vertex] * m_pops[vertex];
            if (accum_rate >= x) {
                delta_pops[out_vertex] = +1;
                return delta_pops;
            }
        }
    }
    return delta_pops;
}

double first_passage_time(gsl_rng *rng, const Model &params, 
        const int &final_vertex) {
    // create an instance of a simulation state:
    gillespie_instance this_run(params);

    while (this_run.m_pops[final_vertex] == 0) {
        this_run.gillespie_step(rng);
    }

    return this_run.m_time;
}
