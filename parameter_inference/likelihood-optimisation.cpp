#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "../graph-model-spec.hpp"
#include "../flying-conjugates.hpp"
#include "../max-likelihood.hpp"
#include "../gillespie-algorithm.hpp"

/* Max. likelihood program:
 *      Generates a simulated dataset with hidden parameter values
 *      Writes out the fake data to a file
 *      Guesses some parameters
 *      Evaluates a likelihood on these parameters
 * Compiles with
 * g++ likelihood-optimisation.cpp ../libs/libgillespie.so ../libs/libflying.so
 *      -lgsl -lm -o guesser
 */

// A higher-order function for mapping a function "mapme" onto data and keeping track of
// the total:
void map_onto_data(Model& params, const epidata_t& this_data, 
                   mappable_t *mapme, real_t *total) {
    // Initial q-values:
    std::vector<real_t> qcorner(params.m_stages, 1); // = {1, 1, 1, 1, 1};
    // integration step
    real_t dt = 0.01;

    for (auto& datum : this_data) {
        // get age, and node id:
        real_t age = datum.first;
        int node = datum.second;
        // and set qvalues accordingly:
        std::vector<real_t> qvals = qcorner;
        // the node we are asking about must be zeroed:
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
        // capture derivative, this is the hazard:
        real_t dprob = prob - prob2;
        dprob /= dt;

        // apply the function to the arguments:
        *total += mapme(age, node, prob, dprob);
    }
}

real_t print_test(real_t age, int node, real_t prob, real_t dprob) {
    // print results:
    std::cout << age << ", " << node << ", ";
    std::cout << prob << ", " << dprob << ", ";
    std::cout << -log(dprob) << std::endl;
    return 1.0;
}

void unit_test(Model& params, const epidata_t& all_times) {
    real_t foo = 0;
    map_onto_data(params, all_times, print_test, &foo);
}

real_t logdprob(real_t age, int node, real_t prob, real_t dprob) {
    // return -log the hazard:
    return -log(dprob);
}

real_t logsurvival(Model& params, int node) {
    // return the survival probability for this node:
    double psurv;
    double alpha = params.m_birth[node];
    double beta = params.m_death[node];
    psurv = alpha / (alpha + beta);
    if (std::isnan(psurv)) psurv = 1.0f;
    return -log(psurv);
}

real_t loglikelihood_test(Model& params, const epidata_t& all_times) {
    real_t energy = 0;
    map_onto_data(params, all_times, logdprob, &energy);
    for (auto& datum : all_times) {
        energy += logsurvival(params, datum.second);
    }
    return energy;
}

std::vector<std::pair<double,int>> generate_dataset(int seed, int runs) {
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
    model.m_birth = {0, 0.05, 0.03, 0, 0};
    model.m_death = {0, 0, 0, 0, 0};
    model.m_initial_pops = {1e2, 0, 0, 0, 0};

    std::vector<int> final_vertices = {3, 4};

    printf("Ground truth:\n");
    printf("  mu = %g\n", mu);
    printf("  rloh = %g\n", rloh);
    printf("  fitness1 = %g\n", model.m_birth[1]);
    printf("  fitness2 = %g\n", model.m_birth[2]);
    printf("  inipop = %g\n", model.m_initial_pops[0]);

    // run some simulations and store the time and final node in
    // all_times:
    std::vector<std::pair<double,int>> all_times;
    times_to_final_vertices(model, seed, runs, final_vertices, all_times);

    return all_times;
}

Model instantiate_model(real_t rloh, real_t mu, real_t fitness1, real_t
fitness2, real_t initialpop) {
    Model params(5);
    params.m_migr[0][1] = mu;
    params.m_migr[0][2] = rloh;
    params.m_migr[1][3] = 0.5 * mu;
    params.m_migr[1][4] = 0.5 * rloh;
    params.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    params.m_birth = {0, fitness1, fitness2, 0, 0};
    params.m_death = {0, 0, 0, 0, 0};
    params.m_initial_pops = {initialpop, 0, 0, 0, 0};

    return params;
}

// SIMULATED ANNEALING STUFF:

double cauchyv(double mean, double width) {
    /* return a Cauchy distributed random variable centred on "mean" with scale
     * coefficient "width" */
    double rand_num = (double)rand() / RAND_MAX;
    double quantile = width * tan(M_PI * (rand_num - 0.5f));
    return mean + quantile;
}

double expv(double scale) {
    /* return an exponentially distributed random variable with mean scale */
    double rand_num = (double)rand() / RAND_MAX;
    return scale * -log(rand_num);
}

Model annealing_min(double (*objective)(Model& model, const epidata_t& dataset),
                    Model initial_guess, const epidata_t& dataset) {
    // Constants:
    const double Tmin = 1e-8; // Minimal temperature
    const double delta = 0.98; // The rate of temperature drop
    const double min_width = 1e-8f;

    srand((unsigned)time(NULL));

    // Initialise variables:
    double Temp = 100;             // Initial temperature 

    Model model = initial_guess;
    Model best_guess = model;
    double best_y = objective(best_guess, dataset);
    double old_best_y = best_y; // Used for adaptive width size

    // Width of neighbourhood:
    double w = Temp; // initial width

    real_t mu   = model.m_migr[0][1];
    real_t rloh = model.m_migr[0][2];
    real_t fitness1 = model.m_birth[1];
    real_t fitness2 = model.m_birth[2];
    real_t inipop   = model.m_initial_pops[0];

    // Simulated annealing process:
    unsigned int iter = 0;
    while (Temp > Tmin) {
        real_t new_mu = expv(mu);
        real_t new_rloh = expv(rloh);
        real_t new_fitness1 = expv(fitness1);
        real_t new_fitness2 = expv(fitness2);
        real_t new_inipop = expv(inipop);
        // DEBUG:
        //std::cout << "Best guesses:" << std::endl;
        //std::cout << "  new_mu = " << new_mu << std::endl;
        //std::cout << "  new_rloh = " << new_rloh << std::endl;
        //std::cout << "  new_fitness1 = " << new_fitness1 << std::endl;
        //std::cout << "  new_fitness2 = " << new_fitness2 << std::endl;
        //std::cout << "  new_inipop = " << new_inipop << std::endl;

        Model new_guess = instantiate_model(new_rloh, new_mu, new_fitness1,
                                new_fitness2, new_inipop);
        
        double y_new = objective(new_guess, dataset);
        if (std::isnan(y_new)) continue;
        double delta_y = y_new - best_y;

        if ((delta_y < 0) || ((rand() / (RAND_MAX + 1.0)) < exp(-delta_y / Temp))) {
            model = new_guess;
            best_guess = model;
            best_y = y_new;
        }
        Temp *= delta;
        /* If best_y has improved by more than the current temperature, reduce
         * the width:*/
        if ((old_best_y - best_y) > Temp) {
            w = 0.5 * (min_width + w); // This will smoothly approach min_width
            old_best_y = best_y;
            printf("Shrinking width\n");
        }
        ++iter;
    }

    printf("System fully cooled after %d iterations\n", iter);
    printf("-log L = %.8g\n", best_y);
    return model;
}

int main(int argc, char* argv[]) {
    // some fake data:
    std::vector<std::pair<double,int>> all_times;
    all_times = generate_dataset(5, 100);

    // save the fake data:
    std::ofstream fakedata;
    fakedata.open("syntheticdata.csv");
    for (auto& entry : all_times) {
        fakedata << entry.first << "," << entry.second << "," << std::endl;
    }
    fakedata.close();

    // Guess some initial model parameters:
    real_t rloh = 0.5e-2;
    real_t mu = 0.5e-3;
    real_t fitness1 = 0.2;
    real_t fitness2 = 0.2;
    real_t initialpop = 100;

    Model guess = instantiate_model(rloh, mu, fitness1, fitness2, initialpop);

    // Try out simulated annealing:
    Model best_guess = annealing_min(loglikelihood_test, guess, all_times);
    mu   = best_guess.m_migr[0][1];
    rloh = best_guess.m_migr[0][2];
    fitness1 = best_guess.m_birth[1];
    fitness2 = best_guess.m_birth[2];
    initialpop   = best_guess.m_initial_pops[0];

    std::cout << "Best guesses:" << std::endl;
    std::cout << "  mu = " << mu << std::endl;
    std::cout << "  rloh = " << rloh << std::endl;
    std::cout << "  fitness1 = " << fitness1 << std::endl;
    std::cout << "  fitness2 = " << fitness2 << std::endl;
    std::cout << "  initialpop = " << initialpop << std::endl;

    return 0;
}
