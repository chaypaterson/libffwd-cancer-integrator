#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"

// A numerical integrator for a four-step branching process model of colorectal
// cancer initiation from Zhang et al. 2022. 
// This integrator is based on integrating a Fourier-space version of the
// Kolmogorov forward equations using the method of characteristics.

std::vector<double> flow(const std::vector<double> &gammas,
                         Model &parameters) {
    // return a vector of flows
    std::vector<double> flux;

    for (size_t vertex = 0; vertex < parameters.m_stages; ++vertex) {
        double rate_of_change = 0;
        rate_of_change += parameters.m_birth[vertex] * (gammas[vertex] - 1) * 
                            gammas[vertex];
        if (parameters.m_death[vertex] > 0) {
            rate_of_change += parameters.m_death[vertex] * 
                                (1.0 / gammas[vertex] - 1) *
                                gammas[vertex];
        }
        for (size_t out_vertex = 0; out_vertex < parameters.m_stages; ++out_vertex) {
            rate_of_change += parameters.m_migr[vertex][out_vertex] * 
                                (gammas[out_vertex] - 1) * gammas[vertex];
        }
        flux.push_back(rate_of_change);
    }

    return flux;
}

void step(std::vector<double> &gammas, double &time, double &dt, 
            Model &parameters) {
    // time step using improved Euler
    std::vector<double> flux = flow(gammas, parameters);
    std::vector<double> gammas2 = gammas;
    for (int vertex = 0; vertex < gammas.size(); ++vertex) {
        gammas2[vertex] += flux[vertex] * dt;
    }
    std::vector<double> flux2 = flow(gammas2, parameters);
    for (int vertex = 0; vertex < gammas.size(); ++vertex) {
        gammas[vertex] += 0.5 * (flux[vertex] + flux2[vertex]) * dt;
    }
    //time += dt;
}

double generating_function(std::vector<double> gammas, 
                           std::vector<double> initial_pops) {
    // compute the value of the generating function at given q values (==gammas)
    // given the initial conditions
    double logpsi = 0;

    if (gammas.size() != initial_pops.size())
        return -1; // an error

    auto q = gammas.begin();
    for (auto& nzero : initial_pops) {
        if (*q > 0)
            logpsi += nzero * log(*q);
        else if (nzero > 0)
            return 0;
        ++q;
    }

    return exp(logpsi);
}

void debug_print(const double &time, const std::vector<double> &qvalues,
                 const Model &parameters) {
    std::cout << time << ", ";
    for (int vertex = 0; vertex < parameters.m_stages; ++vertex) {
        std::cout << qvalues[vertex] << ", ";
    }
    std::cout << std::endl;
}

int main() {
    Model parameters(6);

    // System coefficients:
    parameters.m_migr[0][1] = 2.86e-4;
    parameters.m_migr[1][2] = 1.06e-5;
    parameters.m_migr[2][3] = 9.00e-7;
    parameters.m_migr[3][4] = 1.36e-4;
    parameters.m_migr[4][5] = 4.56e-7;
    parameters.m_birth = {0, 0, 0.2, 0.27, 0.27, 0};
    parameters.m_death = {0, 0, 0, 0, 0, 0};
    parameters.m_initial_pops = {1e8, 0, 0, 0, 0, 0}; 
    
    // we want the probability that site 5 is unoccupied: hence
    std::vector<double> default_qvalues = {1, 1, 1, 1, 1, 1};
    std::vector<std::vector<double>> qvalues;
    for (size_t vertex = 0; vertex < parameters.m_stages; ++vertex) {
        qvalues.push_back(default_qvalues);
        // set this vertices q to zero to get the prob. that this vertex is
        // unoccupied:
        qvalues[vertex][vertex] = 0;
    }
    // TODO we might want to iterate over different sites later
    // list of vectors?

    double time = 0.0;
    const double tmax = 80.0;
    double dt = 1.0;

    while (time < tmax) {
        //debug_print(time, qvalues, parameters);
        // advance the system by smaller, finer steps:
        int subdivision = 20;
        double dt2 = dt / subdivision;
        for (size_t vertex = 0; vertex < parameters.m_stages; ++vertex) {
            for (int i = 0; i < subdivision; ++i) {
                step(qvalues[vertex], time, dt2, parameters);
            } 
        }
        time += dt;

        // TODO iterate over different sites to get the probability
        std::cout << time << ", " ;
        for (size_t vertex = 0; vertex < parameters.m_stages; ++vertex) {
            double prob = generating_function(qvalues[vertex], parameters.m_initial_pops);
            std::cout << 1.0 - prob << ", ";
        }
        std::cout << std::endl;
    }

    //debug_print(time, qvalues, parameters);

    return 0;
}
