#include <vector>
#include <map>
#include <cmath>

#include "graph-model-spec.hpp"
#ifndef FLYING_DEF
#define FLYING_DEF

/* Fast forward integration
 * This is a method for calculating survival and extinction probabilities of
 * birth-death-mutation processes on graphs. This method is based on a
 * generating function representation of the process: this enables the
 * generating function to be computed in terms of the method of characteristics.
 *
 * Initial values are chosen for the conjugate coordinates "q[j]", and these
 * values are evolved along a vector field using Heun's method.
 * 
 * This header file contains the function declarations for the implementations
 * in the fast forward library (fast-forward.cpp).
 */

namespace gmsce {

namespace fast_forward {

// This function computes the right-hand-side of the system of ordinary
// differential equations that govern the characteristics:
std::vector<real_t> rhs_flow(const std::vector<real_t> &qcoords,
                             Model &parameters);

// This function updates the q-coordinates of the characteristic curves using
// Heun's method:
void heun_q_step(std::vector<real_t> &qcoords, const real_t &time, real_t &dt, 
                 Model &parameters);

// Other integration methods:
void improvedeuler_q_step(std::vector<real_t> &qcoords, const real_t &time, real_t &dt, 
                 Model &parameters);
void rungekutta_q_step(std::vector<real_t> &qcoords, const real_t &time, real_t &dt, 
                 Model &parameters);
void ralston_q_step(std::vector<real_t> &qcoords, const real_t &time, real_t &dt, 
                 Model &parameters);
void implicit_q_step(std::vector<real_t> &qcoords, const real_t &time, real_t &dt, 
                 Model &parameters);

// This function computes the value of the generating function at given
// q-coordinates:
real_t generating_function(std::vector<real_t> qcoords, 
                           std::vector<real_t> initial_pops);

}

}

// End of header guard:
#endif
