#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

// A one-hit model with immigration
// -(alpha)-> A -(mu)-> B
//            |
//            V (beta)
//
// A has growth rate beta too.
//
// The point of this experiment is to find out if adjusting beta has an effect
// on the probability of visiting B.

std::vector<int> event(double x, const std::vector<int>& n,
                        const double& Gamma, const double& alpha,
                        const double& beta, const double& mu) {
    // return a set of population changes given x and Gamma
    // TODO the leading tests in the following are redundant because of
    // fall-through
    // TODO can this function be generated programmatically from a graph
    // specification? i.e. at runtime?
    if (x <= alpha) // immigration at fixed rate alpha
        return {+1, 0};

    if ((alpha < x) && (x <= alpha + beta * n[0]))
        return {+1, 0};

    if ((alpha + beta * n[0] < x) && (x <= alpha + 2 * beta * n[0]))
        return {-1, 0};

    // escape at rate mu:
    if ((alpha + 2 * beta * n[0] < x))
        return {-1, +1};

    return {0, 0};
}

double first_passage_time(gsl_rng* rng, const double& alpha,
                        const double& beta, const double& mu) {
    //initialise time
    double t = 0;
    std::vector<int> n = {0, 0};

    // etc.
    double Gamma = 1.0;
    double x = 0.0;
    double Deltat = 0.0;

    while (n[1] == 0) {
        // TODO can this call of event rate be achieved programmatically?
        Gamma = alpha + (2 * beta + mu) * n[0]; 

        // now that we have gamma:
        x = gsl_ran_flat(rng, 0., Gamma);

        // generate an event
        auto dn = event(x, n, Gamma, alpha, beta, mu);
        for (int i = 0; i < n.size(); ++i)
            n[i] += dn[i];

        Deltat = gsl_ran_exponential(rng, 1.0 / Gamma);
        t += Deltat;
    }

    return t;
}

void simulate_runs(int seed, int runs, std::vector<double> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    double alpha = 1.0;
    double beta = 0.3; // adjust me
    double mu = 0.01;

    for (int i = 0; i < runs; ++i) {
        results.push_back(first_passage_time(r, alpha, beta, mu));
    }

    gsl_rng_free(r);
}

int main() {
    // NB this function is completely agnostic
    int num_thr = std::thread::hardware_concurrency();
    std::vector<std::thread> simulations(num_thr);

    int seed = 1;

    std::vector<std::vector<double>> times(num_thr);

    // Run some simulations:
    int runs = 2000;
    for (int i = 0; i < num_thr; ++i) {
        simulations.at(i) = std::thread(simulate_runs, seed + i, runs,
                            std::ref(times[i]));
    }

    for (int i = 0; i < num_thr; ++i) {
        simulations.at(i).join();
    }

    // Print results:
    std::vector<double> all_times;
    for (auto time : times) {
        for (auto t2 : time) {
            all_times.push_back(t2);
        }
    }

    // Sort observed first passage times:
    std::sort(all_times.begin(), all_times.end());

    // build histogram
    for (auto t : all_times) {
        std::cout << t << std::endl;
    }

    std::cout << std::endl;

    return 0;
}
