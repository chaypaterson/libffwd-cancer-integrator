#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

#include "graph-model-spec.hpp"
#include "fast-forward.hpp"

/* A numerical integrator based on the generating function representation of an
 * MVK-type birth-death-mutation process.
 * The model:
 *
 *           (s)
 * 0 -(mu0)-> 1 -(mu1)-> 2
 * i.e.: there are 3 nodes, {0,1,2}. Cells start on node 0 at time t=0, and can
 * migrate from 0 to 1, or 1 to 2, at prescribed rates. Cells on node 1 divide
 * at a rate s.
 *
 * This program calculates the probability that node 2 has been visited using an
 * experimental method based on the method of characteristics.
 */

int main() {
    clonal_expansion::Model two_hit(3);
    // System coefficients:
    two_hit.m_migr[0][1] = 0.001;
    two_hit.m_migr[1][2] = 0.001;
    two_hit.m_birth = {1.0, 1.2, 1.0};
    two_hit.m_death = {1.0, 1.0, 1.0};
    two_hit.m_initial_pops = {1, 0, 0}; 
    
    // We want the probability that site 2 is unoccupied: the corresponding set
    // of q-coordinates is (1,1,0) (see documentation):
    std::vector<double> qvalues = {1, 1, 0};

    double time = 0.0;
    const double tmax = 100.0;
    double dt = 1.0;
    std::cout << time << ", " << 1.0 << ", ";
    std::cout << 0.0 << "," << std::endl;

    while (time < tmax) {
        // Subdivide the time step dt into smaller, finer time steps:
        int subdivision = 16;
        double dt2 = dt / subdivision;
        for (int i = 0; i < subdivision; ++i)
            heun_q_step(qvalues, time, dt2, two_hit);
        // Increment time:
        time += dt;

        double prob = generating_function(qvalues, two_hit.m_initial_pops);
        std::cout << time << ", " << prob << ", ";
        std::cout << 1.0 - prob << "," << std::endl;
    }

    return 0;
}
