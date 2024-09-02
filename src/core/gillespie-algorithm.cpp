#include <cstddef>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <graph-model-spec.hpp>
#include <gillespie-algorithm.hpp>

/* Methods for the gillespie_instance class, which represents a dynamical state
 * of a birth-death process on a graph, as simulated by the Gillespie
 * algorithm.
 */

namespace clonal_expansion {

namespace gillespie_ssa {

// a method to perform one step, changing the populations and advancing time
void gillespie_instance::gillespie_step(gsl_rng *rng) {
    double Gamma = get_gamma();
    double x = gsl_ran_flat(rng, 0, Gamma);
    auto delta_pops = x_to_event(x);

    for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
        // Update populations:
        m_pops[vertex] += delta_pops[vertex];
        // Update m_gamma:
        m_gamma += m_parameters.m_birth[vertex] * delta_pops[vertex];
        m_gamma += m_parameters.m_death[vertex] * delta_pops[vertex];
        for (size_t out_vertex = 0;
                out_vertex < m_vertices; ++out_vertex) {
            m_gamma += m_parameters.m_migr[vertex][out_vertex] * delta_pops[vertex];
        }
    }

    double Delta_t = gsl_ran_exponential(rng, 1.0 / Gamma);

    m_time += Delta_t;
}

// a method to compute Gamma
double gillespie_instance::get_gamma() {
    return m_gamma;
}

// set Gamma to an initial value
void gillespie_instance::set_gamma() {
    m_gamma = 0;
    // for each possible process, increment gamma by the event rate = rate
    // parameter * population on this node
    for (size_t vertex = 0; vertex < m_vertices; ++vertex) {
        m_gamma += m_parameters.m_birth[vertex] * m_pops[vertex];
        m_gamma += m_parameters.m_death[vertex] * m_pops[vertex];
        for (size_t out_vertex = 0;
                out_vertex < m_vertices; ++out_vertex) {
            m_gamma += m_parameters.m_migr[vertex][out_vertex] * m_pops[vertex];
        }
    }
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

// TODO namespace Gillespie {

double first_passage_time(gsl_rng *rng, const Model &params, 
        const int final_vertex) {
    // create an instance of a simulation state:
    gillespie_instance this_run(params);

    // guard against global extinction:
    int total_pop = 0;
    do {
        total_pop = 0;
        for (auto& vertex_pop : this_run.m_pops)
            total_pop += vertex_pop;
        this_run.gillespie_step(rng);
    } while ((this_run.m_pops[final_vertex] == 0) && (total_pop > 0));

    if (total_pop <= 0)
        return -1;

    return this_run.m_time;
}

std::pair<double,int> first_passage_time(gsl_rng *rng, const Model &params,
        const std::vector<int> final_vertices) {
    // create an instance of a simulation state:
    gillespie_instance this_run(params);

    // guard against global extinction:
    int total_pop = 0;
    int ending_type = 0;
    bool repeat = true;
    do {
        total_pop = 0;
        for (auto& vertex_pop : this_run.m_pops)
            total_pop += vertex_pop;
        this_run.gillespie_step(rng);
        for (auto& final_vertex : final_vertices) {
            if (this_run.m_pops[final_vertex] > 0) {
                repeat = false;
                ending_type = final_vertex;
            }
        }
    } while (repeat && (total_pop > 0));

    if (total_pop <= 0)
        return std::make_pair(-1,-1);

    return std::make_pair(this_run.m_time, ending_type);
}

double first_passage_time_poly(gsl_rng *rng, const Model &params,
        const std::vector<int> final_vertices) {
    // create an instance of a simulation state:
    gillespie_instance this_run(params);

    do {
        // guard against global extinction:
        int total_pop = 0;
        for (auto& vertex_pop : this_run.m_pops)
            total_pop += vertex_pop;
        if (!total_pop) break;

        // otherwise, time step:
        this_run.gillespie_step(rng);

        // check if both vertices are occupied
        bool full = 1;
        for (auto& final_vertex : final_vertices) {
            full &= this_run.m_pops[final_vertex];
        }
        if (full) break;
    } while (1);

    return this_run.m_time;
}

void times_to_final_vertex(const Model &model, int seed,
                           int runs_per_thr, int final_vertex,
                           std::vector<double> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    for (int i = 0; i < runs_per_thr; ++i) {
        double time = first_passage_time(r, model, final_vertex);
        if (time >= 0)
            results.push_back(time);
    }

    gsl_rng_free(r);
}

void times_to_final_vertices(const Model &model, int seed,
                             int runs_per_thr, std::vector<int> final_vertices,
                             std::vector<std::pair<double,int>> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    for (int i = 0; i < runs_per_thr; ++i) {
        std::pair<double,int> result =
            first_passage_time(r, model, final_vertices);
        if (result.first >= 0)
            results.push_back(result);
    }
    gsl_rng_free(r);
}

void times_to_final_vertices_poly(const Model &model, int seed,
                             int runs_per_thr, std::vector<int> final_vertices,
                             std::vector<double> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    for (int i = 0; i < runs_per_thr; ++i) {
        double result = first_passage_time_poly(r, model, final_vertices);
        if (result >= 0)
            results.push_back(result);
    }
    gsl_rng_free(r);
}

void print_results(std::vector<double> &all_times) {
    for (auto t : all_times) {
        std::cout << t << std::endl;
    }
}

void print_kaplan_meier(double time_max, std::vector<double> &all_times) {
    size_t num_survivors = all_times.size();
    print_kaplan_meier(time_max, all_times, num_survivors);
}

void print_kaplan_meier(double time_max, std::vector<double> &all_times, size_t ref_pop) {
    size_t num_survivors = ref_pop;
    double dt = time_max / 200;
    double time = 0;
    double survival = 1.0;
    auto time_datum = all_times.begin();
    while (time < time_max) {
        while (*time_datum < time) {
            survival *= (1 - 1.0 / num_survivors);
            --num_survivors;
            ++time_datum;
            // avoid unguarded access:
            if (time_datum >= all_times.end()) break; // break, don't return.
        }
        std::cout << time << ", " << survival << ", ";
        std::cout << 1.0 - survival << "," << std::endl;
        time += dt;
    }
}

double surv_kaplan_meier(double age, std::vector<double> &all_times, size_t ref_pop) {
    // return the value of the Kaplan-Meier estimator at a given age=time
    size_t num_survivors = ref_pop;
    double time_max = all_times.back();
    double dt = time_max / 200;
    double time = 0;
    double survival = 1.0; // TODO this is not very efficient: maybe pass
    // initial survival and age as parameters somehow?
    auto time_datum = all_times.begin();
    while ((time < age) && (time_datum != all_times.end())) {
        // avoid unguarded access:
        if (time >= time_max) break;
        while (*time_datum < time) {
            survival *= (1 - 1.0 / num_survivors);
            --num_survivors;
            ++time_datum;
        }
        time += dt;
    }
    return survival;
}

void print_naive_estimator(double time_max, std::vector<double> &all_times) {
    // ONLY VALID WHEN THERE IS A SINGLE END NODE!
    size_t num_survivors = all_times.size();
    const size_t total_survivors = num_survivors;
    double dt = time_max / 200;
    double time = 0;
    double survival = 1.0;
    auto time_datum = all_times.begin();
    while (time < time_max) {
        while (*time_datum < time) {
            survival = (double)num_survivors / (double)total_survivors;
            --num_survivors;
            ++time_datum;
            // avoid unguarded access:
            if (time_datum == all_times.end()) return;
        }
        std::cout << time << ", " << survival << ", ";
        std::cout << 1.0 - survival << "," << std::endl;
        time += dt;
    }
}

}

}
