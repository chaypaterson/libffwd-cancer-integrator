#include <vector>
#include <map>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

/* Flying conjugates
 * This is a method for calculating survival and extinction probabilities of
 * birth-death-mutation processes on graphs. This method is based on a
 * generating function representation of the process: this enables the
 * generating function to be computed in terms of the method of characteristics.
 *
 * Initial values are chosen for the conjugate coordinates "q[j]", and these
 * values are evolved along a vector field using Heun's method.
 * 
 * This source file contains the function implementations. It has no entry
 * point, and should be compiled to a library with
 *      g++ flying-conjugates.cpp -c -o libflying.o
 */

// TODO more detailed comments!
// TODO maybe too much functionality in each function?
// TODO locate constants centrally

std::vector<double> rhs_flow(const std::vector<double> &qcoords,
                             Model &parameters) {
    // Compute and return a vector containing the rates of change of
    // the characteristic curves (qcoords). This is the right-hand-side of the
    // system of differential equations that govern the characteristics.
    std::vector<double> flux;

    // Sum the different contributions to the total rate of change at each
    // vertex:
    for (size_t vertex = 0; vertex < parameters.m_stages; ++vertex) {
        double rate_of_change = 0;
        // Rate of change of gamma due to birth process:
        rate_of_change += parameters.m_birth[vertex] * 
                            (qcoords[vertex] - 1) * 
                            qcoords[vertex];

        // Rate of change of gamma due to death process:
        rate_of_change += parameters.m_death[vertex] * 
                            (1 - qcoords[vertex]);

        // Rate of change of gamma due to mutation away from this vertex:
        for (size_t out_vertex = 0; 
             out_vertex < parameters.m_stages; 
             ++out_vertex) {
            rate_of_change += parameters.m_migr[vertex][out_vertex] * 
                                (qcoords[out_vertex] - 1) * 
                                qcoords[vertex];
        }

        // Store the total rate of change:
        flux.push_back(rate_of_change);
    }

    return flux;
}

void heun_q_step(std::vector<double> &qcoords, const double &time, double &dt, 
                 Model &parameters) {
    // Time step the q-coordinates (qcoords) using improved Euler integration
    // (also known as Heun's method).
    // Compute an initial guess, qcoords2:
    std::vector<double> flux = rhs_flow(qcoords, parameters);
    std::vector<double> qcoords2 = qcoords;
    for (int vertex = 0; vertex < qcoords.size(); ++vertex) {
        qcoords2[vertex] += flux[vertex] * dt;
    }

    // Use the initial guess qcoords2 to compute an improved guess:
    std::vector<double> flux2 = rhs_flow(qcoords2, parameters);
    // Combine the initial guess and improved guess using the trapezoid rule:
    for (int vertex = 0; vertex < qcoords.size(); ++vertex) {
        qcoords[vertex] += 0.5 * (flux[vertex] + flux2[vertex]) * dt;
    }

    // Do not increment time within this function.
}

double generating_function(std::vector<double> qcoords, 
                           std::vector<double> initial_pops) {
    // Compute the value of the generating function (psi) at given 
    // q-coordinates (qcoords) with known initial conditions.
    // The initial conditions are the initial_pops.

    // To avoid rounding errors when the generating function psi is close to
    // zero or one, first compute the log of psi:
    double log_psi = 0;

    auto q = qcoords.begin();
    for (auto& nzero : initial_pops) {
        if (*q > 0)
            log_psi += nzero * log(*q);
        else if (nzero > 0)
            return 0;
        ++q;
    }

    // Return the exponential of log_psi to get psi:
    return exp(log_psi);
}


