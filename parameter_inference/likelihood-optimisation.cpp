#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

#include "../graph-model-spec.hpp"
#include "../flying-conjugates.hpp"
#include "../max-likelihood.hpp"
#include "../gillespie-algorithm.hpp"

/* Max. likelihood program:
 *      Generates a simulated dataset with hidden parameter values
 *      Writes out the fake data to a file
 *      Guesses some parameters
 *      Evaluates a likelihood on these parameters
 *      Optimises this likelihood with simulated annealing
 * Compiles with
 * g++ likelihood-optimisation.cpp ../libs/libgillespie.so ../libs/libflying.so
 *      -lgsl -lm -o guesser
 */

real_t logsurvival(Model& params, int node) {
    // return the survival probability for this node:
    double psurv;
    double alpha = params.m_birth[node];
    double beta = params.m_death[node];
    psurv = alpha / (alpha + beta);
    if (std::isnan(psurv)) psurv = 1.0f;
    return -log(psurv);
}

real_t loglikelihood_hist_node(Model& params, size_t node, real_t binwidth,
                            size_t ref_population,
                            const std::vector<size_t>& freqs) {
    // Histogram version of the -log likelihood
    // Recieves a histogram of cancers with a known type (node)
    // Returns a -log binomial likelihood
    // TODO fix, pass end_nodes as parameter
    size_t end_nodes[] = {3,4};
    real_t mlogl = 0;

    real_t time = 0;
    std::vector<real_t> qvalsAll(params.m_stages, 1);
    for (auto& endnode : end_nodes)
        qvalsAll[endnode] = 0;
    std::vector<real_t> qvalsExcept = qvalsAll;
    qvalsExcept[node] = 1;
    real_t end_time = binwidth;
    real_t dt = 0.10;
	size_t nsurv = ref_population;

    for (const size_t& curr_bin : freqs) {
        // compute survival probabilities S at start and end of the bin:
        real_t PsiAll = generating_function(qvalsAll, params.m_initial_pops);
        real_t PsiExcept = generating_function(qvalsExcept, params.m_initial_pops);
        real_t Sprob = PsiAll / PsiExcept;
        while (time < end_time) {
            heun_q_step(qvalsAll, time, dt, params);
            heun_q_step(qvalsExcept, time, dt, params);
            time += dt;
        }
        real_t PsiAll2 = generating_function(qvalsAll, params.m_initial_pops);
        real_t PsiExcept2 = generating_function(qvalsExcept, params.m_initial_pops);
        real_t Sprob2 = PsiAll2 / PsiExcept2;

        // -log binomial likelihood:
        real_t p = Sprob - Sprob2;
        mlogl += -log(p) * curr_bin;
        mlogl += -log(1 - p) * (nsurv - curr_bin);
		nsurv -= curr_bin;

        // weight for survival/chance of detection of cancer:
        mlogl += logsurvival(params, node) * curr_bin;
        // update the end time for the next pass:
        end_time += binwidth;
    }

    return mlogl;
}

real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               std::map<size_t, std::vector<size_t>> histos) {
    // Histogram version of the -log likelihood
    // Applied to a set of histograms of cancer incidence with given types
    // The data structure maps the type to a histogram that gives age incidence
    real_t mlogl = 0;

    for (const auto& type : histos) {
        mlogl += loglikelihood_hist_node(params, type.first, binwidth,
                                          ref_population, type.second);
    }
    // TODO end correction: probability not to get either type of cancer
    // some population remaining_pop in the reference pop have not got either
    // cancer. both cancers combined have some lifetime risk lifetime_risk.
    //mlogl += -log(1.0 - lifetime_risk) * remaining_pop;

    return mlogl;
}

std::vector<size_t> convert_to_histogram(const epidata_t& all_times, 
                                         real_t binwidth,
                                         size_t node) {
    real_t max_age = 0;
    for (auto& datum : all_times) {
        if (datum.first > max_age) max_age = datum.first;
    }

    size_t num_bins = max_age / binwidth + 1;

    std::vector<size_t> histogram(num_bins, 0);

    for (auto& entry : all_times) {
        if (entry.second == node) {
            size_t bin = entry.first / binwidth;
            histogram[bin]++;
        }
    }

    return histogram;
}

// Gillespie algorithm functions:

std::vector<std::pair<double,int>> generate_dataset(Model& model, int seed, int runs) {
    std::vector<int> final_vertices = {3, 4};

    // run some simulations and store the time and final node in
    // all_times:
    std::vector<std::pair<double,int>> all_times;
    // TODO namespaces
    times_to_final_vertices(model, seed, runs, final_vertices, all_times);

    return all_times;
}

Model instantiate_model(real_t rloh, real_t mu, real_t fitness1, 
                        real_t fitness2, real_t initialpop) {
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

double logcauchyv(double mode, double width) {
    /* return a log-Cauchy distributed random variable centred on "mode" with scale
     * coefficient "width" */
    double rand_num = (double)rand() / RAND_MAX;
    double quantile = width * tan(M_PI * (rand_num - 0.5f));
    return mode * exp(quantile);
}

double uniform(double mean, double width) {
    double rand_num = (double)rand() / RAND_MAX;
    double deviate = width * 2 * (rand_num - 0.5f);
    return mean + deviate;
}

Model get_neighbour(Model& model, double w) {
    // TODO try pinning some values and fitting others
    real_t new_mu = logcauchyv(model.m_migr[0][1], w);
    real_t new_rloh = logcauchyv(model.m_migr[0][2], w);
    // results are very sensitive to fitness values:
    real_t new_fitness1 = uniform(model.m_birth[1], 0.010 * w);
    real_t new_fitness2 = uniform(model.m_birth[2], 0.010 * w);
    // inipop is unidentifiable:
    real_t new_inipop = model.m_initial_pops[0];

    return instantiate_model(new_rloh, new_mu, new_fitness1,
                            new_fitness2, new_inipop);
}
 
real_t estimate_error(Model& old_model, Model& new_model) {
    // TODO this function is ugly
    // There should be a nice internal method in the Model class for iterating
    // over parameters TODO implement this
    // compute the square % difference in model parameters
    real_t error = 0;
    if (old_model.m_stages != new_model.m_stages) return error;
    for (size_t v_in = 0; v_in < old_model.m_stages; ++v_in) {
        for (size_t v_out = 0; v_out < old_model.m_stages; ++v_out) {
            real_t diff = log(old_model.m_migr[v_in][v_out]) -
                          log(new_model.m_migr[v_in][v_out]);
            if (std::isfinite(diff)) error += diff * diff;
        }

        real_t diff = log(old_model.m_birth[v_in]) -
                      log(new_model.m_birth[v_in]);
        if (std::isfinite(diff)) error += diff * diff;

        diff = log(old_model.m_death[v_in]) -
               log(new_model.m_death[v_in]);
        if (std::isfinite(diff)) error += diff * diff;

        diff = log(old_model.m_initial_pops[v_in]) -
               log(new_model.m_initial_pops[v_in]);
        if (std::isfinite(diff)) error += diff * diff;
    }

    return sqrt(error);
}

/* MAYBE TODO define a class that we can use to "capture" or hide the extraneous
 * arguments like the underlying data and parameters of the histogram. We should
 * be able to pass one object into annealing_min and call it like a function
 * with "objective(model);". This should make the signature of annealing_min
 *      Model annealing_min(ObjectiveFunc objective, Model initial_guess);
 *
 * CURRENTLY I am passing a closure to this object -- this allows the other
 * parameters of the log-likelihood objective function to be "bound" in the
 * context.
 */
Model annealing_min(std::function<real_t(Model& model)> objective, 
                    Model initial_guess) {
    // Function reads a set of datum points and returns a set of model
    // parameters.
    // Constants:
    const double Tmin = 1e-16; // Minimal temperature
    const double delta = 0.98; // The rate of temperature drop
    const double min_width = 1e-2; // min scale of the log-cauchy step
    const double smoothing_factor = 0.82; // for smoothing of width
    const unsigned int iter_max = 2e4;

    // Initialise variables:
    Model best_guess = initial_guess;
    double best_y = objective(best_guess);

    // Width of neighbourhood:
    double w = log(2); // initial width for log-cauchy variates

    printf("Initial likelihood:\n-log L = %g\n", best_y);

    // Simulated annealing process:
    double Temp = best_y;      // Initial temperature 
    unsigned int iter = 0;  // count iterations

    while ((++iter < iter_max) && (Temp > Tmin)) {
        Model new_guess = get_neighbour(best_guess, w);
       
        double y_new = objective(new_guess);

        if (std::isnan(y_new)) continue;

        double delta_y = y_new - best_y;

        if ((delta_y < 0) || ((rand() / (RAND_MAX + 1.0)) < exp(-delta_y / Temp))) {
            // Update best guess:
            best_guess = new_guess;
            best_y = y_new;

            /* Shrink the width: */
            w *= 1.0 - smoothing_factor;
            w += smoothing_factor * min_width;
        }

        Temp *= delta;
    }

    printf("System fully cooled after %d iterations\n", iter);
    printf("-log L = %.8g\n", best_y);
    return best_guess;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Call this program with\n ./guesser seed dataset_size\n");
        return 1;
    }
    // make some fake data:
    std::cout << "Generating synthetic dataset..." << std::endl;
    size_t seed = std::atoi(argv[1]);
    size_t dataset_size = std::atoi(argv[2]);
    srand(seed);

    std::vector<std::pair<double,int>> all_times;
    Model ground_truth = instantiate_model(5.0e-7, 
                                           5.0e-8, 
                                           0.05, 
                                           0.03, 
                                           1e6);

    printf("Ground truth:\n");
    printf("  mu = %g\n", ground_truth.m_migr[0][1]);
    printf("  rloh = %g\n", ground_truth.m_migr[0][2]);
    printf("  fitness1 = %g\n", ground_truth.m_birth[1]);
    printf("  fitness2 = %g\n", ground_truth.m_birth[2]);
    printf("  inipop = %g\n", ground_truth.m_initial_pops[0]);

    all_times = generate_dataset(ground_truth, seed, dataset_size);
    std::cout << "Done. Saving..." << std::endl;

    // compute maximum age:
    real_t max_age = 0;
    {
        // save the fake data:
        std::ofstream fakedata;
        fakedata.open("syntheticdata_raw.csv");
        for (auto& entry : all_times) {
            fakedata << entry.first << "," << entry.second << "," << std::endl;
            max_age += (entry.first > max_age) * (entry.first - max_age);
        }
        fakedata.close();
    }

    // Convert age data to histogram:
    size_t reference_pop = all_times.size();
    real_t binwidth = max_age / (2 * pow(reference_pop, 0.4)); // years
    std::map<size_t, std::vector<size_t>> incidence;
    incidence[3] = convert_to_histogram(all_times, binwidth, 3);
    incidence[4] = convert_to_histogram(all_times, binwidth, 4);
    printf("Target likelihood:\n-log L = %g\n", 
            loglikelihood_hist_both(ground_truth, binwidth, reference_pop,
                                    incidence));

    // Save the histogram:
    {
        std::ofstream fakedata;
        fakedata.open("syntheticdata_hist.csv");
        fakedata << "bin width = " << binwidth << "\n";
        fakedata << "max age = " << max_age << "\n";
        fakedata << "ref population = " << reference_pop << std::endl;

        fakedata << "[";
        for (auto& bin : incidence[3]) 
            fakedata << bin << ", ";
        fakedata << "]" << std::endl;
        fakedata << "[";
        for (auto& bin : incidence[4])
            fakedata << bin << ", ";
        fakedata << "]" << std::endl;

        fakedata.close();
    }

    unsigned tries = 0, maxtries = 10;
    {
        // Guess some initial model parameters:
        real_t rloh = 1e-7;
        real_t mu = 1e-8;
        real_t fitness1 = 0.05;
        real_t fitness2 = 0.03;
        real_t initialpop = 1e6;

        Model guess = instantiate_model(rloh, mu, fitness1, fitness2, initialpop);

        // Try out simulated annealing:
        std::cout << "Starting annealing..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth, 
                                           reference_pop, incidence);
            // the histogram version
        };

        Model best_guess = annealing_min(objective, guess);
        // Annealing now complete. Print guess:
        std::cout << "Best guesses:" << std::endl;
        std::cout << "  mu = " << best_guess.m_migr[0][1] << std::endl;
        std::cout << "  rloh = " << best_guess.m_migr[0][2] << std::endl;
        std::cout << "  fitness1 = " << best_guess.m_birth[1] << std::endl;
        std::cout << "  fitness2 = " << best_guess.m_birth[2] << std::endl;
        std::cout << "  initialpop = " << best_guess.m_initial_pops[0] << std::endl;
    }

    return 0;
}
