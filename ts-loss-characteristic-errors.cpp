#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

/* A conjugate characteristics simulation of tumour suppressor loss
 * This program estimates the global error with Richardson extrapolation.
 */

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "please provide a time step: \n";
        std::cout << "  ./this_program time_step";
        std::cout << std::endl;
        return 1;
    }

    // System coefficients:
    real_t rloh = 0.5e-2;
    real_t mu = 0.5e-3;

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
    std::vector<real_t> qvalues = {1, 1, 1, 1, 0};
    std::vector<real_t> qvalues2 = qvalues;

    real_t time = 0.0;
    const real_t tmax = 100.0;
    real_t dt = atof(argv[1]); // integration step
    real_t half = 0.5 * dt;
    real_t t_write_step = 1.0; // write out step

    // print key:
    std::cout << "age, err," << std::endl;
    std::cout.precision(20);
    // print first line:
    std::cout << std::fixed << time << ", ";
    std::cout << std::fixed << 0.0 << ", " << std::endl;
    real_t t_write = time + t_write_step;

    while (time < tmax) {
        while (time < t_write - dt) {
            heun_q_step(qvalues, time, dt, model);
            // Advance qvalues 2 by half the relevant time step, twice
            heun_q_step(qvalues2, time, half, model);
            heun_q_step(qvalues2, time, half, model);
            // Increment time:
            time += dt;
        }

        {
            real_t delta = t_write - time;
            heun_q_step(qvalues, time, delta, model);
            real_t halfdelta = 0.5 * delta;
            heun_q_step(qvalues2, time, halfdelta, model);
            heun_q_step(qvalues2, time, halfdelta, model);
            time = t_write;
            // compute the P(t) values:
            real_t prob = generating_function(qvalues, model.m_initial_pops);
            real_t prob2 = generating_function(qvalues2, model.m_initial_pops);
            // compute the corresponding error:
            double err = (prob - prob2);
            err /= (2 << 2 - 1); // richardson extrapolation
            err *= (1 - 2 * (err < 0)); // abs value
            // Write out errors:
            std::cout << std::fixed << time << ", " << err << std::endl;
            t_write += t_write_step;
        }
    }

    return 0;
}
