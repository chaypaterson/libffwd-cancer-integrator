#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>

// A numerical integrator based on the generating function representation of an
// MVK-type birth-death-mutation process
// The model:
//
//           (s)
// A -(mu0)-> B -(mu1)-> C

std::vector<double> flow(const std::vector<double> &gammas,
                        const double &mu0, 
                        const double &s, const double &mu1) {
    // return a vector of flows
    std::vector<double> flux = {0,0,0};

    flux[0] = mu0 * (gammas[1] - 1) * gammas[0];
    flux[1] = s * (gammas[1] - 1) * gammas[1];
    flux[1] += mu1 * (gammas[2] - 1) * gammas[1];

    return flux;
}

void step(std::vector<double> &gammas, double &time, double &dt,
        const double &mu0, const double &s, const double &mu1) {
    // time step using improved Euler
    std::vector<double> flux = flow(gammas, mu0, s, mu1);
    std::vector<double> gammas2 = gammas;
    for (int i = 0; i < gammas.size(); ++i) {
        gammas2[i] += flux[i] * dt;
    }
    std::vector<double> flux2 = flow(gammas2, mu0, s, mu1);
    for (int i = 0; i < gammas.size(); ++i) {
        gammas[i] += 0.5 * (flux[i] + flux2[i]) * dt;
    }

    time += dt;
}

double generating_function(std::vector<double> gammas, 
                           std::vector<int> initial_pops) {
    // compute the value of the generating function at given q values (==gammas)
    // given the initial conditions
    double logpsi = 0;

    if (gammas.size() != initial_pops.size())
        return -1; // an error

    auto q = gammas.begin();
    for (auto& nzero : initial_pops) {
        if (*q > 0)
            logpsi += nzero * log(*q);
        ++q;
    }

    return exp(logpsi);
}

int main() {
    // System coefficients:
    const double mu0 = 0.001;
    const double mu1 = 0.001;
    const double s = 0.2;
    
    // we want the probability that site 2 is unoccupied: hence
    std::vector<double> qvalues = {1, 1, 0};

    std::vector<int> initial_pops = {1, 0, 0}; 

    double time = 0.0;
    const double tmax = 100.0;
    double dt = 1.0;
    std::cout << time << ", " << 1.0 << ", ";
    std::cout << 0.0 << "," << std::endl;

    while (time < tmax) {
        // advance the system by smaller, finer steps:
        int subdivision = 20;
        double dt2 = dt / subdivision;
        for (int i = 0; i < subdivision; ++i)
            step(qvalues, time, dt2, mu0, s, mu1);
        // beautiful

        double prob = generating_function(qvalues, initial_pops);
        std::cout << time << ", " << prob << ", ";
        std::cout << 1.0 - prob << "," << std::endl;
    }

    return 0;
}
