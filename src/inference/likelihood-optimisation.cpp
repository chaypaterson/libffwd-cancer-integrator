#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "graph-model-spec.hpp"
#include "fast-forward.hpp"
#include "gillespie-algorithm.hpp"
#include <thread>

/* Max. likelihood program:
 *      Generates a simulated dataset with hidden parameter values
 *      Writes out the fake data to a file
 *      Guesses some parameters
 *      Evaluates a likelihood on these parameters
 *      Optimises this likelihood with simulated annealing
 * Compiles with
 * g++ likelihood-optimisation.cpp ../libs/libgillespie.so ../libs/libflying.so
 *      -lgsl -lm -o guesser
 *
 * Contact:
 *      Chay Paterson (chay.paterson@manchester.ac.uk)
 *      Miaomiao Gao (miaomiao.gao@postgrad.manchester.ac.uk)
 */

using clonal_expansion::gillespie_ssa::times_to_final_vertices;
using clonal_expansion::fast_forward::generating_function;
using clonal_expansion::fast_forward::heun_q_step;
using clonal_expansion::real_t;
using clonal_expansion::Model;

typedef std::vector<std::pair<real_t, int>> epidata_t;

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

real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               std::map<size_t, std::vector<size_t>> histos_sporadic,
                               std::map<size_t, std::vector<size_t>> histos_germline) {
    // this version includes two separate sporadic and germline incidence curves.
    real_t mlogl = 0;

    for (const auto& type : histos_sporadic) {
        mlogl += loglikelihood_hist_node(params, type.first, binwidth,
                                         ref_population, type.second);
    }

    Model params_germline = params;
    params_germline.m_initial_pops[1] = params.m_initial_pops[0];
    params_germline.m_initial_pops[0] = 0;

    for (const auto& type : histos_germline) {
        mlogl += loglikelihood_hist_node(params, type.first, binwidth,
                                         ref_population, type.second);
    }

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
    times_to_final_vertices(model, seed, runs, final_vertices, all_times);
    // TODO this is a natural candidate for parallelisation

    return all_times;
}

Model instantiate_model(real_t rloh, real_t mu, real_t fitness1,
                        real_t fitness2, real_t initialpop) {
    // Spawn a model of tumour suppressor loss
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

Model instantiate_model_germline(real_t rloh, real_t mu, real_t fitness1,
                                 real_t fitness2, real_t initialpop) {
    // Spawn a model of tumour suppressor loss
    Model params(5);
    params.m_migr[0][1] = mu;
    params.m_migr[0][2] = rloh;
    params.m_migr[1][3] = 0.5 * mu;
    params.m_migr[1][4] = 0.5 * rloh;
    params.m_migr[2][4] = 0.5 * mu;
    // birth and death rates:
    params.m_birth = {0, fitness1, fitness2, 0, 0};
    // TODO maybe fitness1 should be 0 for this model? discuss w colleagues
    params.m_death = {0, 0, 0, 0, 0};
    params.m_initial_pops = {0, initialpop, 0, 0, 0};
    // TODO germline mutations: which initial node? i think 1?

    return params;
}

Model differ_model(Model& params, real_t drloh, real_t dmu, real_t dfitness1,
                   real_t dfitness2, real_t dinitialpop) {
    Model dmodel = params;
    dmodel.m_migr[0][1] += dmu;
    dmodel.m_migr[0][2] += drloh;
    dmodel.m_migr[1][3] += 0.5 * dmu;
    dmodel.m_migr[1][4] += 0.5 * drloh;
    dmodel.m_migr[2][4] += 0.5 * dmu;
    // birth and death rates:
    dmodel.m_birth[1] += dfitness1;
    dmodel.m_birth[2] += dfitness2;
    dmodel.m_initial_pops[0] += dinitialpop;

    return dmodel;
}

Model differ_model(Model& params, std::vector<real_t> Delta) {
    return differ_model(params, Delta[0], Delta[1], Delta[2], Delta[3], 0);
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
    // results are very sensitive to fitness values (so use a narrower distribution):
    real_t new_fitness1 = uniform(model.m_birth[1], 0.10 * w);
    real_t new_fitness2 = uniform(model.m_birth[2], 0.10 * w);
    // force the fitnesses to be positive (negative values are meaningless in this context)
    new_fitness1 *= 1 - 2 * (new_fitness1 < 0);
    new_fitness2 *= 1 - 2 * (new_fitness2 < 0);
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
Model annealing_min(std::function<real_t(Model &model)> objective,
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

// Numerical analysis methods:

// Stencil for numerical differentiation:
std::vector<std::vector<double>> StencilLP16() {
    /* forms a 16 point stencil that combines a derivative with a low pass
     * filter using a tensor product
     */
    std::vector<std::vector<double>> stencil;

    std::vector<double> coeffs = {-0.125, -0.25, 0, 0.25, 0.125};
    for (int m = -2; m <= 2; ++m) {
        for (int n = -2; n <= 2; ++n) {
            double weight = coeffs[m + 2] * coeffs[n + 2];
            if (weight * weight > 0)
                stencil.push_back({weight, (double)m, (double)n});
        }
    }

    return stencil;
}

Eigen::MatrixXd compute_gradient(std::function<real_t(Model&)> objective,
                                 Model point) {
    // Compute the gradient of the objective function at a point with respect to
    // the logs of the model parameters:
    int dim = 4;
    Eigen::MatrixXd Gradient(dim,1);

    // Numerical derivatives:
    double epsilon = 1e-3; // use the same epsilon as we used for the hessian

    std::vector<double> weights = {-0.125, -0.25, 0, 0.25, 0.125};
    int radius = (weights.size() - 1) / 2; // NB: weights.size must be odd

    // vector of parameter values:
    std::vector<real_t> Theta = {point.m_migr[0][2] /*rloh*/,
                                 point.m_migr[0][1] /*mu*/,
                                 point.m_birth[1]   /*fitness1*/,
                                 point.m_birth[2]   /*fitness2*/
                                };

    for (int axis = 0; axis < dim; ++axis) {
        /* Compute the derivative along each axis using a finite difference
           scheme and the weights given above: */
        double diff = 0;
        for (int tap = 0, step = -radius; tap < weights.size(); ++tap, ++step) {
            std::vector<real_t> Delta(dim, 0);
            Delta[axis] = epsilon * Theta[axis] * step;
            Model dmodel = differ_model(point, Delta);
            diff += weights[tap] * objective(dmodel) / epsilon;
        }

        Gradient(axis, 0) = diff;
    }

    return Gradient;
}

Eigen::MatrixXd compute_hessian(std::function<real_t(Model&)> objective,
                                Model point) {
    // Compute the hessian of the objective function at a point:
    int dim = 4;
    Eigen::MatrixXd Hessian(dim,dim);

    // Numerical derivatives:
    double epsilon = 4e-3;
    /* Try to choose epsilon to balance truncation error and
     * catastrophic cancellations.
     *
     * - Chay
     */
    // vector of parameter values:
    std::vector<real_t> Theta = {point.m_migr[0][2] /*rloh*/,
                                 point.m_migr[0][1] /*mu*/,
                                 point.m_birth[1]   /*fitness1*/,
                                 point.m_birth[2]   /*fitness2*/
                                };

    // Pre-compute a stencil and weights to use for finite differencing:
    std::vector<std::vector<double>> stencil;
    stencil = StencilLP16();

    // To improve numerical stability, take the derivatives with regards to
    // the logs of the parameters, then convert back:
    for (int x = 0; x < dim; ++x) {
        for (int y = 0; y < dim; ++y) {
            // Compute the Hessian using a finite difference scheme:

            double diff = 0;
            for (auto& tap : stencil) {
                std::vector<real_t> Delta(dim, 0);
                Delta[x] += epsilon * Theta[x] * tap[1];
                Delta[y] += epsilon * Theta[y] * tap[2];
                Model dmodel = differ_model(point, Delta);
                diff += tap[0] * objective(dmodel) / epsilon;
            }
            // Scale diff appropriately to get second derivative:
            diff /= epsilon;

            Hessian(x,y) = diff;
        }
    }

    // Raw Hessian (in nice units):
    std::cout << "H = " << std::endl;
    std::cout << "[rloh, mu, s1, s2]" << std::endl;
    std::cout << Hessian << std::endl;

    // Convert from log coords to true Hessian:
    Eigen::MatrixXd Jacobian(dim,dim);
    for (int i = 0; i < dim; ++i)
        Jacobian(i,i) = 1.0 / Theta[i];

    Hessian = Jacobian.transpose() * Hessian * Jacobian;

    return Hessian;
}

Model gradient_min(std::function<real_t(Model& model)> objective,
                   Model initial_guess) {
    /* Minimise objective using gradient descent (NB: we are using numerical
       differentiation, not autodiff/backpropagation)
     */

    // Initialise variables:
    Model best_guess = initial_guess;
    double best_y = objective(best_guess);

    int dim = 4;
    double learning_rate = 1e-4;
    double tolerance = 1e-3; // % change tolerated in parameter estimates

    // while the gradient is too large:
    while (1) {
        // compute the gradient at the current best guess
        Eigen::MatrixXd gradient = compute_gradient(objective, best_guess);

        std::vector<real_t> Theta = {best_guess.m_migr[0][2] /*rloh*/,
                                     best_guess.m_migr[0][1] /*mu*/,
                                     best_guess.m_birth[1]   /*fitness1*/,
                                     best_guess.m_birth[2]   /*fitness2*/
                                    };

        // update the current best guess by -learning_rate * gradient
        std::vector<real_t> Delta(dim, 0);
        for (int axis = 0; axis < dim; ++axis) {
            Delta[axis] = Theta[axis] * (exp(-learning_rate * gradient(axis)) - 1);
            std::cout << Delta[axis] << std::endl;
        }

        best_guess = differ_model(best_guess, Delta);

        // we want the percentage-change in all coefficients to be lower than
        // tolerance:

        Eigen::MatrixXd Jacobian(dim,dim);
        double change = 0;
        for (int i = 0; i < dim; ++i) {
            Jacobian(i,i) = 1.0 / Theta[i];
            double dx = Delta[i] / Theta[i];
            change += dx * dx;
        }

        change = std::sqrt(change);

        // TODO DEBUG:
        std::cout << change << std::endl;

        // if the gradient is small enough, exit:
        if (change < tolerance) break;
        if (std::isnan(change)) break;
    }

    return best_guess;
}

// Methods to output data:

void print_model(Model &model) {
    printf("  mu = %g\n", model.m_migr[0][1]);
    printf("  rloh = %g\n", model.m_migr[0][2]);
    printf("  fitness1 = %g\n", model.m_birth[1]);
    printf("  fitness2 = %g\n", model.m_birth[2]);
    printf("  inipop = %g\n", model.m_initial_pops[0]);
}

void write_model_line(std::ofstream& file, Model &model) {
    file << model.m_migr[0][1];
    file << ", ";
    file << model.m_migr[0][2];
    file << ", ";
    file << model.m_birth[1];
    file << ", ";
    file << model.m_birth[2];
    file << ", ";
    file << model.m_initial_pops[0];
    file << ",";
    file << std::endl;
}

void save_histogram(real_t& binwidth, real_t&  max_age, size_t& reference_pop,
                    std::vector<size_t>& end_nodes,
                    std::map<size_t, std::vector<size_t>>& incidence,
                    std::string filename) {
    std::ofstream fakedata;
    fakedata.open(filename);
    fakedata << "bin width = " << binwidth << "\n";
    fakedata << "max age = " << max_age << "\n";
    fakedata << "ref population = " << reference_pop << std::endl;

    for (auto &end_node : end_nodes) {
        fakedata << "[";
        for (auto& bin : incidence[end_node])
            fakedata << bin << ", ";
        fakedata << "]" << std::endl;
    }

    fakedata.close();
}

void save_histogram(real_t& binwidth, real_t&  max_age, size_t& reference_pop,
                    std::vector<size_t>& end_nodes,
                    std::map<size_t, std::vector<size_t>>& incidence) {
    save_histogram(binwidth, max_age, reference_pop, end_nodes, incidence,
                   "syntheticdata_hist.csv");
}

std::map<size_t, std::vector<size_t>> jackknife_incidence(
                                       size_t index,
                                       const std::map<size_t, std::vector<size_t>>
                                       &histogram,
                                       std::vector<size_t> end_nodes) {
    // delete the index-th entry in the histogram
    std::map<size_t, std::vector<size_t>> new_incidence = histogram;
    // first we have to find the indexth person
    size_t count = 0;
    for (auto& end_node : end_nodes) {
        for (auto& bin : new_incidence[end_node]) {
            count += bin;
            if (count > index) {
                // delete one individual:
                --bin;
                // we are now done:
                return new_incidence;
            }
        }
    }
    // if we get here something has gone wrong are there are not enough entries
    // in the histogram.
    return new_incidence;
}

void resample_incidence(
    const std::map<size_t, std::vector<size_t>> *incidence,
    size_t reference_pop, size_t start, size_t end, std::vector<size_t> end_nodes,
    real_t binwidth, Model *initial_guess, std::vector<Model> *resampled_estimates) {

    for (unsigned tries = start; tries < end; ++tries) {
        // Resample the incidence:
        std::map<size_t, std::vector<size_t>> resampled_incidence;
        resampled_incidence = jackknife_incidence(tries, *incidence, end_nodes);

        std::cout << "Resampling..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, resampled_incidence);
            // the histogram version
        };

        Model best_guess = annealing_min(objective, *initial_guess);
        // Annealing now complete. Print guess:
        std::cout << "Best guesses:" << std::endl;
        print_model(best_guess);
        // Annealing now complete. Store guessed model parameters:
        resampled_estimates->push_back(best_guess);
    }
}

void jackknife_and_save(std::map<size_t, std::vector<size_t>> &incidence,
                        size_t reference_pop, real_t binwidth,
                        std::vector<size_t> end_nodes, Model &initial_guess,
                        size_t n_child_threads = 0) {
    // spawn child threads to carry out resampling
    // how many resamples should we generate on each child thread?
    // in total we need reference_pop many, but simply running
    // reference_pop/(n_child_threads+1) replicates on each thread might run too few
    // if reference_pop is not a round multiple of n_child_threads. instead, run
    // reference_pop/(n_child_threads+1) on all child threads, and
    // the remainder, reference_pop % n_child_threads, in the
    // main thread. this should also allow this program to run in
    // single-threaded mode in the same way, making it backwards-compatible
    size_t runs_per_thr = reference_pop / (n_child_threads + 1);
    size_t remainder    = reference_pop % (n_child_threads + 1);

    std::vector<std::vector<Model>> results(n_child_threads + 1);

    std::vector<std::thread> child_threads(0);
    for (int thread = 0; thread < n_child_threads; ++thread) {
        size_t start, end; // start and end replicates: range of data-points in
                           // reference population to resample
        start = thread * runs_per_thr;
        end = start + runs_per_thr;

        child_threads.push_back(std::thread(resample_incidence,
            &incidence, reference_pop, start, end, end_nodes, binwidth,
            &initial_guess, &results[thread]));
    }

    // run runs_per_thr + remainder in this, the parent thread:
    resample_incidence(&incidence, reference_pop, reference_pop - runs_per_thr - remainder,
                       reference_pop, end_nodes, binwidth,
                       &initial_guess, &results[n_child_threads]);

    // Wait for child threads to finish
    for (auto& thread : child_threads) {
        thread.join();
    }

    // save the point estimates so we can make a density plot later
    std::ofstream estimates_by_row;
    estimates_by_row.open("resampled_estimates.csv");
    estimates_by_row << "mu, rloh, s1, s2, initial_pop," << std::endl;
    for (auto& estimate_set : results) {
        for (auto& estimate : estimate_set) {
            write_model_line(estimates_by_row, estimate);
        }
        estimate_set.clear();
    }
    estimates_by_row.close();
}


int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Call this program with\n ./guesser seed dataset_size\n");
        return 1;
    }
    size_t seed = std::atoi(argv[1]);
    size_t dataset_size = std::atoi(argv[2]);
    bool include_germline = 0; // run with germline or without?
    srand(seed);

    std::vector<std::pair<double,int>> all_times, all_times_germline;

    Model ground_truth = instantiate_model(5.0e-7,
                                           5.0e-8,
                                           0.05,
                                           0.03,
                                           1e6);

    printf("Ground truth:\n");
    print_model(ground_truth);

    // simulate some data that could be collected by a longitudinal study:
    std::cout << "Generating synthetic dataset..." << std::endl;
    if (include_germline) {
        // if we include germline cases, the clincal study has a 50:50 mix of
        // sporadic cases and cases with germline alterations
        Model ground_truth_germline = instantiate_model_germline(5.0e-7,
                                      5.0e-8,
                                      0.05,
                                      0.03,
                                      1e6);
        printf("Ground truth (germline):\n");
        print_model(ground_truth_germline);

        all_times = generate_dataset(ground_truth, seed, dataset_size / 2);
        // to include germline mutations:
        // vary the initial populations so that instead of node 0 the cells are
        // initially on node 1
        all_times_germline = generate_dataset(ground_truth_germline, seed,
                                              dataset_size / 2);
    } else {
        // sporadic cases only (default)
        all_times = generate_dataset(ground_truth, seed, dataset_size);
    }

    std::vector<size_t> end_nodes = {3,4};
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
    size_t reference_pop = all_times.size(); /* NB: with germline this is 1/2
        the previous value */
    real_t binwidth = max_age / (2 * pow(reference_pop, 0.4)); // years

    std::map<size_t, std::vector<size_t>> incidence, incidence_germline;

    for (auto &end_node : end_nodes) {
        incidence[end_node] = convert_to_histogram(all_times, binwidth,
                              end_node);
        if (include_germline) {
            incidence_germline[end_node] = convert_to_histogram(
                                               all_times_germline, binwidth,
                                               end_node);
        }
    }

    // Save the histogram:
    save_histogram(binwidth, max_age, reference_pop, end_nodes, incidence);
    if (include_germline) {
        save_histogram(binwidth, max_age, reference_pop, end_nodes,
                       incidence_germline, "syntheticdata_hist_germline.csv");
    }

    // Get the main un-resampled estimate with simulated annealing:
    Model best_guess = ground_truth;
    {
        std::cout << "--------------" << std::endl;
        std::cout << "Best estimate:" << std::endl;
        // Guess some initial model parameters:
        real_t rloh = 1e-7;
        real_t mu = 1e-8;
        real_t fitness1 = 0.05;
        real_t fitness2 = 0.03;
        real_t initialpop = 1e6;

        if (include_germline) {
            printf("Target likelihood:\n-log L = %g\n",
                   loglikelihood_hist_both(ground_truth, binwidth, reference_pop,
                                           incidence, incidence_germline));
        } else {
            printf("Target likelihood:\n-log L = %g\n",
                   loglikelihood_hist_both(ground_truth, binwidth, reference_pop,
                                           incidence));
        }

        Model guess = instantiate_model(rloh, mu, fitness1, fitness2, initialpop);

        // Get and return Hessian:
        Eigen::MatrixXd Hessian;

        // Try out simulated annealing:
        std::cout << "Starting annealing..." << std::endl;
        std::function<real_t(Model&)> objective;
        if (include_germline) {
            objective = [&](Model& model) {
                return loglikelihood_hist_both(model, binwidth,
                                               reference_pop, incidence,
                                               incidence_germline);
            };
            // this should now work with germline data
        } else {
            objective = [&](Model& model) {
                return loglikelihood_hist_both(model, binwidth,
                                               reference_pop, incidence);
            };
        }
        best_guess = annealing_min(objective, guess);
        Hessian = compute_hessian(objective, best_guess);

        // Annealing now complete. Print guess:
        std::cout << "Best guesses:" << std::endl;
        print_model(best_guess);


        std::cout << "H = " << std::endl;
        std::cout << "[rloh, mu, s1, s2]" << std::endl;
        std::cout << Hessian << std::endl;

        printf("\n");

        // Invert Hessian to get covariance matrix:
        std::cout << "cov = H^-1 = " << std::endl;
        std::cout << Hessian.inverse() << std::endl;

        // Compute confidence intervals:
        std::cout << "Standard deviations:" << std::endl;
        int param = 0;
        for (int param = 0; param < 4; ++param) {
            double estimate;
            // TODO disgusting:
            if (param == 0) estimate = best_guess.m_migr[0][2];
            if (param == 1) estimate = best_guess.m_migr[0][1];
            if (param == 2) estimate = best_guess.m_birth[1];
            if (param == 3) estimate = best_guess.m_birth[2];
            std::cout << estimate << " +/- ";
            std::cout << sqrt(Hessian.inverse()(param, param)) << std::endl;
        }
    }

    // Jack-knife resampling of incidence, now using gradient descent instead of
    // simulated annealing:
    jackknife_and_save(incidence, reference_pop, binwidth, end_nodes, best_guess);

    return 0;
}
