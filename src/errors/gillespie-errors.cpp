#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

#include <gillespie-algorithm.hpp>

// A Gillespie algorithm simulation of tumour suppressor loss
// This program returns expected errors in the Gillespie simulation

using clonal_expansion::gillespie_ssa::surv_kaplan_meier; 
using clonal_expansion::gillespie_ssa::times_to_final_vertices; 
using clonal_expansion::real_t;

std::map<int,std::vector<double>> generate_dataset(int seed, int runs) {
    // System coefficients:
    double rloh = 5.0e-7;
    double mu = 5.0e-8;

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
        return 1;
    }

    int seed = atoi(argv[1]);
    int runs = atoi(argv[2]);

    std::map<int,std::vector<double>> all_times_1 = generate_dataset(seed, runs);
    std::map<int,std::vector<double>> all_times_2 = generate_dataset(seed + runs, runs);

    // compute maximum age in either dataset, and sort both datasets:
    double age_max = 0;
    for (auto& pair : all_times_1) {
        for (auto& age : pair.second) {
            age_max += (age > age_max) * (age - age_max);
        }
        std::sort((pair.second).begin(), (pair.second).end());
    }

    for (auto& pair : all_times_2) {
        for (auto& age : pair.second) {
            age_max += (age > age_max) * (age - age_max);
        }
        std::sort((pair.second).begin(), (pair.second).end());
    }

    // compute scaled root mean square difference in survival curves:
    size_t reference_pop_1 = all_times_1[3].size() + all_times_1[4].size();
    size_t reference_pop_2 = all_times_2[3].size() + all_times_2[4].size();
    // NB: these should both be equal to runs

    size_t num_sample_points = 256;
    double dt = age_max / (double)num_sample_points;

    for (double age = 0; age <= age_max; age += dt) {
        double toterr = 0;
        for (int type = 3; type < 5; ++type) {
            real_t s1 = 1, s2 = 1;
            if (all_times_1[type].size())
                s1 = surv_kaplan_meier(age, all_times_1[type], reference_pop_1);
            if (all_times_2[type].size())
                s2 = surv_kaplan_meier(age, all_times_2[type], reference_pop_2);
            real_t error = s1 - s2;
            error *= error;
            toterr += error;
        }
        toterr /= 2;
        toterr = std::sqrt(toterr);
        std::cout.precision(9);
        std::cout << std::fixed << age << ",";
        std::cout.precision(20);
        std::cout << std::fixed << toterr << "," << std::endl;
    }

    return 0;
}
