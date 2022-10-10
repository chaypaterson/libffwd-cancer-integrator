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

void simulate_runs(const Model &model, int seed, 
            int runs_per_thr, int final_vertex,
            std::vector<double> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    for (int i = 0; i < runs_per_thr; ++i) {
        results.push_back(first_passage_time(r, model, final_vertex));
    }

    gsl_rng_free(r);
}

void print_results(std::vector<double> &all_times) {
    for (auto t : all_times) {
        std::cout << t << std::endl;
    }
}

void print_kaplan_meier(double time_max, std::vector<double> &all_times) {
    size_t num_survivors = all_times.size();
    double dt = time_max / 100;
    double time = 0;
    double survival = 1.0;
    auto time_datum = all_times.begin();
    while (time < time_max) {
        while (*time_datum < time) {
            survival *= (1 - 1.0 / num_survivors);
            --num_survivors;
            ++time_datum;
        }
        std::cout << time << ", " << survival << ", ";
        std::cout << 1.0 - survival << "," << std::endl;
        time += dt;
    }
}

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
    std::vector<double> all_times;

    // run (num_thr * runs_per_thr) simulations and store the times in
    // all_times:
    {
        std::vector<std::vector<double>> times(num_thr);
        {
            std::vector<std::thread> simulations(num_thr);

            // Run some simulations:
            for (int i = 0; i < num_thr; ++i) {
                simulations.at(i) = std::thread(simulate_runs, 
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
