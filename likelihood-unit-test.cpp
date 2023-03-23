#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"

/* Max. likelihood unit test
 * Test definition of likelihood function
 */

void unit_test(Model& params, std::vector<std::pair<real_t, int>>& all_times) {
    // Initial q-values:
    std::vector<real_t> qcorner(params.m_stages, 1); // = {1, 1, 1, 1, 1};
    // integration step
    real_t dt = 0.01;

    for (auto& datum : all_times) {
        // get age, and node id:
        real_t age = datum.first;
        int node = datum.second;
        // and set qvalues accordingly:
        std::vector<real_t> qvals = qcorner;
        qvals[node] = 0;
        // integrate to get likelihood:
        real_t time = 0.0;
        while (time < age) {
            heun_q_step(qvals, time, dt, params);
            time += dt;
        }

        real_t prob = generating_function(qvals, params.m_initial_pops);
        // advance one dt step into the future:
        heun_q_step(qvals, time, dt, params);
        real_t prob2 = generating_function(qvals, params.m_initial_pops);
        // capture derivative:
        real_t dprob = prob - prob2;
        dprob /= dt;

        // print results:
        std::cout << age << ", " << node << ", ";
        std::cout << prob << ", " << dprob << ", ";
        std::cout << -log(dprob) << std::endl;
    }
}

real_t loglikelihood_test(Model& params, 
                          std::vector<std::pair<real_t, int>>& all_times) {
    // Initial q-values:
    std::vector<real_t> qcorner(params.m_stages, 1); // = {1, 1, 1, 1, 1};
    // integration step
    real_t dt = 0.01;
    real_t energy = 0;

    for (auto& datum : all_times) {
        // get age, and node id:
        real_t age = datum.first;
        int node = datum.second;
        // and set qvalues accordingly:
        std::vector<real_t> qvals = qcorner;
        qvals[node] = 0;
        // integrate to get likelihood:
        real_t time = 0.0;
        while (time < age) {
            heun_q_step(qvals, time, dt, params);
            time += dt;
        }

        real_t prob = generating_function(qvals, params.m_initial_pops);
        // advance one dt step into the future:
        heun_q_step(qvals, time, dt, params);
        real_t prob2 = generating_function(qvals, params.m_initial_pops);
        // capture derivative:
        real_t dprob = prob - prob2;
        dprob /= dt;

        energy += -log(dprob);
    }

    return energy;
}

int main(int argc, char* argv[]) {
    // System coefficients:
    real_t rloh = 0.5e-2;
    real_t mu = 0.5e-3;

    Model params(5);
    params.m_migr[0][1] = mu;
    params.m_migr[0][2] = rloh;
    params.m_migr[1][3] = 0.5 * mu;
    params.m_migr[1][4] = 0.5 * rloh;
    params.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    params.m_birth = {0, 0.2, 0.2, 0, 0};
    params.m_death = {0, 0, 0, 0, 0};
    params.m_initial_pops = {1e2, 0, 0, 0, 0};

    std::vector<std::pair<double,int>> all_times;

    // some fake data:
    all_times.push_back(std::make_pair(33.0, 3));
    all_times.push_back(std::make_pair(50.0, 4));

    // pass params, num. nodes and all_times into unit test:
    unit_test(params, all_times);

    // second unit test:
    std::cout << loglikelihood_test(params, all_times);
    std::cout << std::endl;

    return 0;
}
