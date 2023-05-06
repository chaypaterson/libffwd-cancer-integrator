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

std::vector<std::pair<double,int>> generate_dataset(int seed, int runs) {
    // System coefficients:
    double rloh = 0.5e-2;
    double mu = 0.5e-3;

    Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0.05, 0.03, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e2, 0, 0, 0, 0};

    std::vector<int> final_vertices = {3, 4};

    // run some simulations and store the time and final node in
    // all_times:
    std::vector<std::pair<double,int>> all_times;
    times_to_final_vertices(model, seed, runs, final_vertices, all_times);

    return all_times;
}

int main(int argc, char* argv[]) {
    int seed = atoi(argv[1]);
    int runs = atoi(argv[2]);

    std::vector<std::pair<double,int>> all_times = generate_dataset(seed, runs);

    std::cout << "age, node," << std::endl;
    for (auto& pair : all_times) {
        std::cout << pair.first << ", " << pair.second << "," << std::endl;
    }

    return 0;
}
