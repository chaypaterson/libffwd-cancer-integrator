#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

// A conjugate characteristics simulation of tumour suppressor loss

int main(int argc, char* argv[]) {
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
    // TODO new conditions: what about the probability to end up on 3, the
    // probability to end up on 3 or 4, and the probability to end up on 4 GIVEN
    // THAT 3 or 4 have been visited?
    std::vector<real_t> qvalues = {1, 1, 1, 1, 0};

    real_t time = 0.0;
    const real_t tmax = 100.0;
    real_t dt = atof(argv[1]); // integration step
    real_t t_write_step = 1.0; // write out step

    // print key:
    std::cout << "age, p1, p2," << std::endl;
    std::cout.precision(17);
    // print first line:
    std::cout << std::fixed << time << ", " << 1.0 << ", ";
    std::cout << std::fixed << 0.0 << ", " << std::endl;
    real_t t_write = time + t_write_step;

    while (time < tmax) {
        heun_q_step(qvalues, time, dt, model);
        // Increment time:
        time += dt;

        real_t prob = generating_function(qvalues, model.m_initial_pops);
        if (time >= t_write) {
            // If we have overshot the write out time, go back a bit:
            real_t delta = t_write - time;
            heun_q_step(qvalues, time, delta, model);
            time = t_write;
            // Write out probabilities:
            std::cout << std::fixed << time << ", " << prob << ", ";
            std::cout << std::fixed << 1.0 - prob << ", " << std::endl;
            t_write += t_write_step;
        }
    }

    return 0;
}
