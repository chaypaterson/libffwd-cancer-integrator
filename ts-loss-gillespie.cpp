#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

#include "gillespie-algorithm.hpp"

// A Gillespie algorithm simulation of tumour suppressor loss

int main() {
    int num_thr = std::thread::hardware_concurrency();
    int runs_per_thr = 1e7;
    int seed = 19;

    // System coefficients:
    double rloh = 1e-4;
    double mu = 1e-6;

    Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0, 0, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e3, 0, 0, 0, 0};

    std::vector<int> final_vertices = {3, 4};

    // Run many simulations and store the results:
    std::vector<std::pair<double,int>> all_times;

    // run (num_thr * runs_per_thr) simulations and store the times in
    // all_times:
    {
        std::vector<std::vector<std::pair<double,int>>> times(num_thr);
        {
            std::vector<std::thread> simulations(num_thr);

            // Run some simulations:
            for (int i = 0; i < num_thr; ++i) {
                // Use a different seed for each simulation:
                simulations.at(i) = std::thread(times_to_final_vertices, 
                                    model, seed + i, runs_per_thr,
                                    final_vertices, 
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

    std::vector<double> mutant_times;
    for (auto& result : all_times)
        if (result.second == 4)
            mutant_times.push_back(result.first);
    // Kaplan-Meier plot:
    std::cout << "age, p1, p2," << std::endl;
    print_kaplan_meier(100, mutant_times);

    std::cout << std::endl;

    return 0;
}
