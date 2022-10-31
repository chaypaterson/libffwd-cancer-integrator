#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

// A conjugate characteristics simulation of tumour suppressor loss

std::vector<real_t> energy_coeffs(const Model& model) {
    // compute multipliers for the energy computation:
    // The multipliers should be zero on nodes with no birth, non-zero
    // everywhere else, and sum to zero:
    std::vector<real_t> coefficients(model.m_stages,0);
    int num_nonzero_nodes = 0;
    for (int vertex = 0; vertex < model.m_stages; ++vertex) {
        if (model.m_birth[vertex] > 0) {
            num_nonzero_nodes++;
            coefficients[vertex] = -1;
        }
    }
    for (int vertex = 0; vertex < model.m_stages; ++vertex) {
        if (model.m_birth[vertex] > 0) {
            coefficients[vertex] += num_nonzero_nodes;
            break;
        }
    }

    return coefficients;
}

real_t energy_guess(std::vector<real_t> qvalues, const Model& model) {
    // compute a conserved quantity for this model and q-trajectory
    real_t energy = 0;
    std::vector<real_t> coefficients = energy_coeffs(model);
    // compute an initial guess:
    for (int vertex = 0; vertex < model.m_stages; ++vertex) {
        if (model.m_birth[vertex] > 0) {
            energy += coefficients[vertex] * 
                (log(1 - qvalues[vertex]) - 
                    log(qvalues[vertex] * model.m_birth[vertex] 
                        - model.m_death[vertex]))
                / (model.m_birth[vertex] - model.m_death[vertex]);
        }
    }

    return energy;
}

real_t energy_update(real_t& energy, std::vector<real_t> qvalues, real_t dt, Model& model) {
    // updates an "energy" checksum. In theory this function should evaluate to
    // zero, so the accumulated value in energy gives an indicator of the
    // accumulated numerical errors.
    std::vector<real_t> coefficients = energy_coeffs(model);
    // compute the correction:
    real_t correction = 0;
    std::vector<real_t> rate_of_change = rhs_flow(qvalues, model);
    for (int vertex = 0; vertex < model.m_stages; ++vertex) {
        if (model.m_birth[vertex] > 0) {
            real_t denominator = (model.m_birth[vertex] * qvalues[vertex] -
                                    model.m_death[vertex]) * 
                                 (qvalues[vertex] - 1);
            correction += coefficients[vertex] * dt *
                          rate_of_change[vertex] / denominator;
            for (int out_vtx = 0; out_vtx < model.m_stages; ++out_vtx) {
                correction -= coefficients[vertex] * dt * 
                              model.m_migr[vertex][out_vtx] *
                              (qvalues[out_vtx] - 1) * qvalues[vertex] 
                                / denominator;
            }
        }
    }
    energy += correction;
    return correction;
}

int main() {
    // System coefficients:
    real_t rloh = 0.5e-2;
    real_t mu = 0.5e-3;

    Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0.2, 0.2, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e2, 0, 0, 0, 0};

    // final vertices 3 and 4: 4 = mutants with LOH
    std::vector<real_t> qvalues = {1, 1, 1, 1, 0};

    real_t time = 0.0;
    const real_t tmax = 100.0;
    real_t dt = 0.5;
    // initialize energy to zero (arbitrary value):
    real_t energy = 0.0;
    // print key:
    std::cout << "age, p1, p2, energy, delta," << std::endl;
    // print first line:
    std::cout << time << ", " << 1.0 << ", ";
    std::cout << 0.0 << ", " << energy << ", " ;
    std::cout << 0.0 << "," << std::endl;

    while (time < tmax) {
        // Subdivide the time step dt into smaller, finer time steps:
        int subdivision = 16;
        real_t dt2 = dt / subdivision;
        for (int i = 0; i < subdivision; ++i)
            heun_q_step(qvalues, time, dt2, model);
        // Increment time:
        time += dt;
        // energy drift correction:
        if (std::isnan(energy)) energy = energy_guess(qvalues, model);
        real_t delta = energy_update(energy, qvalues, dt, model);

        real_t prob = generating_function(qvalues, model.m_initial_pops);
        std::cout << time << ", " << prob << ", ";
        std::cout << 1.0 - prob << ", " << energy << ", " << delta;
        std::cout << "," << std::endl;
    }

    return 0;
}
