#include <cstdio>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <thread>
#include <algorithm>

// A Gillespie algorithm simulation of a birth-death-mutation process
// Also a two-hit model of cancer initiation
// The model:
//
//           (s)
// A -(mu0)-> B -(mu1)-> C
//
// compile me with g++ two-hit.cpp -lgsl -pthread

std::vector<int> event(double x, const std::vector<int> &n, 
                        const double &Gamma, const double &mu0, 
                        const double &s, const double &mu1) {
    // return a set of population changes
    if (x <= mu0 * n[0])
        return {0, +1, 0};
    if ((mu0 * n[0] < x) && (x <= (mu0 * n[0] + s * n[1])))
        return {0, +1, 0};
    if ((mu0 * n[0] + s * n[1]) < x)
        return {0, 0, +1};

    return {0,0,0};
}

double first_passage_time(gsl_rng *rng, const double &mu0, const double &s,
    const double &mu1) {
    // initialise time
    double t = 0.;
    // populations
    std::vector<int> n = {1, 0, 0};

    // etc.
    double Gamma = 1.0;
    double x = 0.0;
    double Deltat = 0.0;

    // main loop
    while (n[2] == 0) {
        Gamma = mu0 * n[0] + (mu1 + s) * n[1];
        x = gsl_ran_flat(rng, 0., Gamma);

        // get event
        auto dn = event(x, n, Gamma, mu0, s, mu1);
        for (int i = 0; i < n.size(); ++i)
            n[i] += dn[i];

        Deltat = gsl_ran_exponential(rng, 1.0 / Gamma);
        t += Deltat;
    }

    // give us this passage time
    return t;
}

void simulate_runs(int seed, int runs_per_thr, std::vector<double> &results) {
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // Seed RNG:
    gsl_rng_set(r,seed);

    // System coefficients:
    double mu0 = 0.001;
    double mu1 = 0.001;
    double s = 0.1;

    for (int i = 0; i < runs_per_thr; ++i) {
        results.push_back(first_passage_time(r, mu0, s, mu1));
    }

    gsl_rng_free(r);
}

void print_results(std::vector<double> &all_times) {
    for (auto t : all_times) {
        std::cout << t << std::endl;
    }
}

void print_kaplan_meier(double time_max, std::vector<double> &all_times) {
    size_t num_survivors = all_times.size();
    double dt = time_max / 100;
    double time = 0;
    double survival = 1.0;
    auto time_datum = all_times.begin();
    while (time < time_max) {
        while (*time_datum < time) {
            survival *= (1 - 1.0 / num_survivors);
            --num_survivors;
            ++time_datum;
        }
        std::cout << time << ", " << survival << ", ";
        std::cout << 1.0 - survival << "," << std::endl;
        time += dt;
    }
}

int main() {
    int num_thr = std::thread::hardware_concurrency();
    int runs_per_thr = 1e6;
    int seed = 1;

    // I hypothesise that the timescale for the late anomaly should be around
    // ~50 years. Beyond this point, the distribution should be dominated by
    // mu0, as this is the rate limiting process.
    std::vector<double> all_times;

    // run (num_thr * runs_per_thr) simulations and store the times in
    // all_times:
    {
        std::vector<std::vector<double>> times(num_thr);
        {
            std::vector<std::thread> simulations(num_thr);

            // Run some simulations:
            for (int i = 0; i < num_thr; ++i) {
                simulations.at(i) = std::thread(simulate_runs, seed + i,
                                    runs_per_thr, std::ref(times[i]));
            }

            for (int i = 0; i < num_thr; ++i) {
                simulations.at(i).join();
            }
        }

        // Flatten and store results:
        for (auto time : times) {
            for (auto t2 : time) {
                all_times.push_back(t2);
            }
        }
    }

    std::sort(all_times.begin(), all_times.end());

    // Kaplan-Meier plot:
    print_kaplan_meier(100, all_times);

    std::cout << std::endl;

    return 0;
}
