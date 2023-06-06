#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

#include "gillespie-algorithm.hpp"

// A Gillespie algorithm simulation of a birth-death-mutation process
// Also a two-hit model of cancer initiation
// The model:
//
//           (s)
// A -(mu0)-> B -(mu1)-> C
//
// compile me with g++ two-hit.cpp -lgsl -pthread

int main() {
    int num_thr = std::thread::hardware_concurrency();
    int runs_per_thr = 1e7;
    num_thr = 1;
    int seed = 1;

    // System coefficients:
    // these should be packaged in params or Model or
    // something similar.
    double mu0 = 0.001;
    double mu1 = 0.001;

    Model model(3);
    model.m_migr[0][1] = mu0;
    model.m_migr[1][2] = mu1;
    model.m_birth = {1, 1.2, 1};
    model.m_death = {1, 1.0, 1};
    model.m_initial_pops = {1, 0, 0};

    // I hypothesise that the timescale for the late anomaly should be around
    // ~50 years. Beyond this point, the distribution should be dominated by
    // mu0, as this is the rate limiting process.
    // TODO can stash this in a function as it's agnostic to the structure of
    // model:
    std::vector<double> all_times;

    // run (num_thr * runs_per_thr) simulations and store the times in
    // all_times:
    {
        std::vector<std::vector<double>> times(num_thr);
        {
            std::vector<std::thread> simulations(num_thr);

            // Run some simulations:
            for (int i = 0; i < num_thr; ++i) {
                // Use a different seed for each simulation:
                simulations.at(i) = std::thread(times_to_final_vertex, 
                                    model, seed + i, runs_per_thr, 2,
                                    std::ref(times[i]));
            }

            for (int i = 0; i < num_thr; ++i) {
                simulations.at(i).join();
            }
        }

        // Flatten and store results:
        for (auto time : times) {
            for (auto t2 : time) {
                all_times.push_back(t2);
            }
        }
    }

    std::sort(all_times.begin(), all_times.end());

    // Kaplan-Meier plot:
    print_kaplan_meier(100, all_times);

    std::cout << std::endl;

    return 0;
}
