#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

#include "../gillespie-algorithm.hpp"

// A Gillespie algorithm simulation of tumour suppressor loss
// This will generate data for the maximum likelihood section

std::map<int,std::vector<double>> generate_dataset(int seed, int runs) {
    // System coefficients:
    double rloh = 5.0e-7;
    double mu = 5.0e-8;

    Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0.05, 0.03, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e6, 0, 0, 0, 0};

    std::vector<int> final_vertices = {3, 4};

    // run some simulations and store the time and final node in
    // all_times:
    std::vector<std::pair<double,int>> all_times;
    times_to_final_vertices(model, seed, runs, final_vertices, all_times);

    std::map<int,std::vector<double>> all_times_flipped;
    for (auto& entry : all_times) {
        all_times_flipped[entry.second].push_back(entry.first);
    }

    return all_times_flipped;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Call this program with: " << std::endl;
        std::cout << "    ./gillespie-sampler seed runs";
        std::cout << std::endl;
    }

    int seed = atoi(argv[1]);
    int runs = atoi(argv[2]);

    std::map<int,std::vector<double>> all_times = generate_dataset(seed, runs);

    std::cout << "age, node," << std::endl;
    double age_max = 0;
    for (auto& pair : all_times) {
        for (auto& age : pair.second) {
            std::cout << age << ", " << pair.first << "," << std::endl;
            age_max += (age > age_max) * (age - age_max);
        }
        std::sort((pair.second).begin(), (pair.second).end());
    }
    std::cout << std::endl;

    std::cout << "\nSurvival curves:" << std::endl;
    size_t reference_pop = all_times[3].size() + all_times[4].size();

    for (auto& pair : all_times) {
        std::cout << "Type: " << pair.first << std::endl;
        print_kaplan_meier(age_max, pair.second, reference_pop);
    }

    return 0;
}
