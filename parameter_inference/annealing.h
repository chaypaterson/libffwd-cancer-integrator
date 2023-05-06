#ifndef ANNEAL_DEF
#define ANNEAL_DEF

#include <math.h>

// This is the objective function I want to minimise:
double energy(double x);

// This function will accept an objective function, an initial guess for the
// value, and a target precision. It will return the best guess for the minimum
// of obj_fun, to within a precision of "precision".
// I am leaving the exact implementation up to you!
double annealing_min(double (*obj_fun)(double), 
                     double x_initial_guess, 
                     double precision);

// This function describes the probability that we accept a new neighbour.
double prob_accept(double x_new, double x_best, double temperature);

// This is the main function.
// I am leaving the exact implementation up to you!
int main(int argc, char* argv[]);

#endif
