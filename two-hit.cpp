#include <cstdio>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// compile me with g++ two-hit.cpp -lgsl

double first_passage_time(int seed, double mu0, double s, double mu1) {
    // TODO pass the RNG in as an argument?
    // Initialize RNG
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // initialise time
    double t = 0.;
    // populations
    int n0 = 1;
    int n1 = 0;
    int n2 = 0;

    // etc.
    double Gamma = 1.0;
    double x = 0.0;
    double Deltat = 0.0;

    // main loop
    while (n2 == 0) {
        Gamma = mu0 * n0 + (mu1 + s) * n1;
        x = gsl_ran_flat(r, 0., Gamma);
        // get event
        if (x <= mu0 * n0) {
            n0 -= 1;
            n1 += 1;
        }
        if ((mu0 * n0 < x) && (x <= (mu0 * n0 + s * n1))) {
            n1 += 1;
        }
        if ((mu0 * n0 + s * n1) < x) {
            n1 -= 1;
            n2 += 1;
        }

        Deltat = gsl_ran_exponential(r, 1.0 / Gamma);
        t += Deltat;
    }

    gsl_rng_free(r);
    // give us this passage time
    return t;
}

int main() {
    int seed = 1;
    // System coefficients:
    double mu0 = 0.001;
    double mu1 = 0.001;
    double s = 0.1;
    // I hypothesise that the timescale for the late anomaly should be around
    // 50 years.

    // Run some simulations:
    double t2 = first_passage_time(seed, mu0, s, mu1);
    // Print results:
    std::cout << t2 << std::endl;

    return 0;
}
