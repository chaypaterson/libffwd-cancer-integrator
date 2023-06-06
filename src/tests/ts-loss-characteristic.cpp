#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "fast-forward.hpp"

// A conjugate characteristics simulation of tumour suppressor loss

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("provide dt\n");
        return 1;
    }
    // System coefficients:
    real_t rloh = 6.957e-7;
    real_t mu = 5.9606e-8;

    Model model(5);
    model.m_migr[0][1] = mu;
    model.m_migr[0][2] = rloh;
    model.m_migr[1][3] = 0.5 * mu;
    model.m_migr[1][4] = 0.5 * rloh;
    model.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    model.m_birth = {0, 0.04766, 0.03484, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e6, 0, 0, 0, 0};

    // final vertices 3 and 4: 4 = mutants with LOH
    // new conditions: now trying to compute the probability of PRIMARY cases of
    // 3 (i.e. survival for 3 CONDITIONED on 4 being empty)
    // for one set of qvalues, all end nodes should be zero. for the other set,
    // all end nodes should be zero EXCEPT the node of interest.
    std::vector<real_t> qvaluesF = {1, 1, 1, 0, 0};
    std::vector<real_t> qvalues3 = {1, 1, 1, 0, 1};

    real_t time = 0.0;
    const real_t tmax = 380.0;
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
        while (time < t_write - dt) {
            heun_q_step(qvaluesF, time, dt, model);
            heun_q_step(qvalues3, time, dt, model);
            // Increment time:
            time += dt;
        }

        {
            real_t delta = t_write - time;
            heun_q_step(qvaluesF, time, delta, model);
            heun_q_step(qvalues3, time, delta, model);
            time = t_write;
            // Write out probabilities:
            real_t probF = generating_function(qvaluesF, model.m_initial_pops);
            real_t prob3 = generating_function(qvalues3, model.m_initial_pops);
            real_t prob = probF / prob3;
            std::cout << std::fixed << time << ", " << prob << ", ";
            std::cout << std::fixed << 1.0 - prob << ", " << std::endl;
            t_write += t_write_step;
        }
    }

    return 0;
}
