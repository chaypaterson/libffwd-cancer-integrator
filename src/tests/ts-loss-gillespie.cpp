#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

#include <gillespie-algorithm.hpp>

// A Gillespie algorithm simulation of tumour suppressor loss

using clonal_expansion::gillespie_ssa::times_to_final_vertices;
using clonal_expansion::gillespie_ssa::print_kaplan_meier;

int main(int argc, char* argv[]) {
    int runs_per_thr = 1e7; // default values
    int seed = 1;
    if (argc < 3) {
        printf("Call this program with\n ./tsgillespie seed runs \n");
        //printf("type should be 3 or 4\n");
        return 1;
    } // else:
    int num_thr = std::thread::hardware_concurrency() - 2;
    seed = atoi(argv[1]);
    runs_per_thr = atoi(argv[2]) / num_thr;
    //int type = atoi(argv[3]);

    // System coefficients:
    double rloh = 5e-7;
    double mu = 5e-8;

    clonal_expansion::Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0.05, 0.03, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e6, 0, 0, 0, 0};

    // We want the chance that the first tumour detected is of type "type".
    // We want to stop the simulations when the first tumour is detected: i.e.
    // of any type. And then we want the fraction of these that have the chosen
    // type.
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
                simulations.at(i) = std::thread(
                                    times_to_final_vertices, 
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
    size_t study_population = all_times.size();
    double time_max = all_times.back().first;

    // Print results for both types:
    for (auto type : final_vertices) {
        std::cout << "Type " << type << ":" << std::endl;

        std::vector<double> mutant_times;
        for (auto& result : all_times) {
            if (result.second == type) {
                mutant_times.push_back(result.first);
            }
        }

        // Guard against invalid access:
        if (mutant_times.size() < 1) {
            std::cout << "No results" << std::endl;
            return 2;
        }

        // Kaplan-Meier plot:
        std::cout << "age, p1, p2," << std::endl;
        print_kaplan_meier(time_max, mutant_times, study_population);

        std::cout << std::endl;
    }

    return 0;
}
