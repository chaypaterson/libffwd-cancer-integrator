#include <vector>
#include <cmath>

/* Compute a conserved quantity for the flying-conjugates method. This can be
 * used as a checksum in simulations to study the accumulation of numerical
 * errors.
 */

#include "flying-conjugates.hpp"

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
