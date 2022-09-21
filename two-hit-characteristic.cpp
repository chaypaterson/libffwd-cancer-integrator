#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"

// A numerical integrator based on the generating function representation of an
// MVK-type birth-death-mutation process
// The model:
//
//           (s)
// A -(mu0)-> B -(mu1)-> C

// TODO separate out the relevant functions into their own .cpp library!

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
    time += dt;
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

int main() {
    Model two_hit(3);
    // System coefficients:
    two_hit.m_migr[0][1] = 0.001;
    two_hit.m_migr[1][2] = 0.001;
    two_hit.m_birth = {0, 0.2, 0};
    two_hit.m_death = {0, 0, 0};
    two_hit.m_initial_pops = {1, 0, 0}; 
    
    // we want the probability that site 2 is unoccupied: hence
    std::vector<double> qvalues = {1, 1, 0};

    double time = 0.0;
    const double tmax = 100.0;
    double dt = 1.0;
    std::cout << time << ", " << 1.0 << ", ";
    std::cout << 0.0 << "," << std::endl;

    while (time < tmax) {
        // advance the system by smaller, finer steps:
        int subdivision = 20;
        double dt2 = dt / subdivision;
        for (int i = 0; i < subdivision; ++i)
            step(qvalues, time, dt2, two_hit);
        // beautiful

        double prob = generating_function(qvalues, two_hit.m_initial_pops);
        std::cout << time << ", " << prob << ", ";
        std::cout << 1.0 - prob << "," << std::endl;
    }

    return 0;
}
