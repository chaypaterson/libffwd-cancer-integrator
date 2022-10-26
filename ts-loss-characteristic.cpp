#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

// A conjugate characteristics simulation of tumour suppressor loss

int main() {
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
    model.m_birth = {0, 0.2, 0.2, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e2, 0, 0, 0, 0};

    // final vertices 3 and 4: 4 = mutants with LOH
    std::vector<double> qvalues = {1, 1, 1, 1, 0};

    double time = 0.0;
    const double tmax = 100.0;
    double dt = 1.0;
    std::cout << "age, p1, p2," << std::endl;
    std::cout << time << ", " << 1.0 << ", ";
    std::cout << 0.0 << "," << std::endl;

    while (time < tmax) {
        // Subdivide the time step dt into smaller, finer time steps:
        int subdivision = 16;
        double dt2 = dt / subdivision;
        for (int i = 0; i < subdivision; ++i)
            heun_q_step(qvalues, time, dt2, model);
        // Increment time:
        time += dt;

        double prob = generating_function(qvalues, model.m_initial_pops);
        std::cout << time << ", " << prob << ", ";
        std::cout << 1.0 - prob << "," << std::endl;
    }

    return 0;
}
