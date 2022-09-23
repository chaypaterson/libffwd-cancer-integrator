#include <vector>
#include <map>
#include <cmath>

#include "graph-model-spec.hpp"
#ifndef FLYING_DEF
#define FLYING_DEF

/* Flying conjugates
 * This is a method for calculating survival and extinction probabilities of
 * birth-death-mutation processes on graphs. This method is based on a
 * generating function representation of the process: this enables the
 * generating function to be computed in terms of the method of characteristics.
 *
 * Initial values are chosen for the conjugate coordinates "q[j]", and these
 * values are evolved along a vector field using Heun's method.
 * 
 * This header file contains the function declarations for the implementations
 * in the flying-conjugates library (flying-conjugates.cpp).
 */

// This function computes the right-hand-side of the system of ordinary
// differential equations that govern the characteristics:
std::vector<double> rhs_flow(const std::vector<double> &qcoords,
                             Model &parameters);

// This function updates the q-coordinates of the characteristic curves using
// Heun's method:
void heun_q_step(std::vector<double> &qcoords, const double &time, double &dt, 
                 Model &parameters);

// This function computes the value of the generating function at given
// q-coordinates:
double generating_function(std::vector<double> qcoords, 
                           std::vector<double> initial_pops);

// End of header guard:
#endif
