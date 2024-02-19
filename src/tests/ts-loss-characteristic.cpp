#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "fast-forward.hpp"

// A conjugate characteristics simulation of tumour suppressor loss

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printf("Call this program with\n ./tsconj dt type\n");
        printf("type should be 3 or 4\n");
        return 1;
    } // else:
    // System coefficients:
    real_t rloh = 5e-7;
    real_t mu = 5e-8;

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

    // final vertices 3 and 4: 4 = mutants with LOH, 3 = mutants without.
    // new conditions: now trying to compute the probability of PRIMARY cases of
    // 3 (i.e. survival for 3 CONDITIONED on 4 being empty)
    // for one set of qvalues, all end nodes should be zero. for the other set,
    // all end nodes should be zero EXCEPT the node of interest.

    std::vector<real_t> qvaluesBoth = {1, 1, 1, 0, 0};
    std::vector<real_t> qvaluesOther = qvaluesBoth;
    int type = atoi(argv[2]); // type of tumour we are interested in
    qvaluesOther[type] = 1; // prob. that there are no tumours of the other type

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
            fast_forward::heun_q_step(qvaluesBoth, time, dt, model);
            fast_forward::heun_q_step(qvaluesOther, time, dt, model);
            // Increment time:
            time += dt;
        }

        {
            real_t delta = t_write - time;
            fast_forward::heun_q_step(qvaluesBoth, time, delta, model);
            fast_forward::heun_q_step(qvaluesOther, time, delta, model);
            time = t_write;
            // Write out probabilities:
            real_t probneither = fast_forward::generating_function(
                                    qvaluesBoth, 
                                    model.m_initial_pops);
            real_t probother = fast_forward::generating_function(
                                    qvaluesOther, 
                                    model.m_initial_pops);
            // S(type, age) = Pr(neither type, age) / Pr(none of other type,
            // age)
            real_t prob = probneither / probother;
            std::cout << std::fixed << time << ", " << prob << ", ";
            std::cout << std::fixed << 1.0 - prob << ", " << std::endl;
            t_write += t_write_step;
        }
    }

    return 0;
}
