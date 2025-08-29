#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include <thread>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <graph-model-spec.hpp>
#include <fast-forward.hpp>
#include <gillespie-algorithm.hpp>

#include "guesser-argparse.hpp"

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
 *
 * TODO this file is way, way too long and needs to be broken into smaller
 * units. Also, eigen dependencies are very spread out through this file.
 */

using clonal_expansion::gillespie_ssa::times_to_final_vertices;
using clonal_expansion::fast_forward::generating_function;
using clonal_expansion::fast_forward::heun_q_step;
using clonal_expansion::real_t;
using clonal_expansion::Model;
using clonal_expansion::GuesserConfig;

typedef std::vector<std::pair<double, int>> Epidata_t;
typedef std::map<size_t, std::vector<size_t>> Histogram_t;

// Statistical functions:

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
        real_t p = (Sprob - Sprob2) / Sprob; // thanks Luo Xiang Ge!
        if (p > 0 && p < 1) {
            mlogl += -log(p) * curr_bin;
            mlogl += -log(1 - p) * (nsurv - curr_bin);
        }
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
                               Histogram_t histos) {
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
                               Histogram_t histos_sporadic,
                               Histogram_t histos_germline) {
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

std::vector<size_t> convert_to_histogram(const Epidata_t& all_times,
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

Epidata_t generate_dataset(Model& model, int seed, int runs) {
    std::vector<int> final_vertices = {3, 4};

    // run some simulations and store the time and final node in
    // all_times:
    Epidata_t all_times;
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

std::vector<real_t> model_params_pure(Model point) {
    // return a vector of values that summarise the model parameters
    std::vector<real_t> params = {point.m_migr[0][2] /*rloh*/,
                                  point.m_migr[0][1] /*mu*/,
                                  point.m_birth[1]   /*fitness1*/,
                                  point.m_birth[2]   /*fitness2*/
                                 };
    return params;
}

std::vector<volatile real_t*> model_params_raw(Model point) {
    // return a vector of raw pointers to the model parameters
    std::vector<volatile real_t*> params = {&point.m_migr[0][2] /*rloh*/,
                                            &point.m_migr[0][1] /*mu*/,
                                            &point.m_birth[1]   /*fitness1*/,
                                            &point.m_birth[2]   /*fitness2*/
                                           };
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
    // TODO maybe fitness1 should be 0 for this model? discuss w Miriam and
    // David
    params.m_death = {0, 0, 0, 0, 0};
    params.m_initial_pops = {0, initialpop, 0, 0, 0};

    return params;
}

Model shifted_model(Model& params, real_t drloh, real_t dmu, real_t dfitness1,
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

Model shifted_model(Model& params, std::vector<real_t> Delta) {
    return shifted_model(params, Delta[0], Delta[1], Delta[2], Delta[3], 0);
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
    real_t new_fitness1 = logcauchyv(model.m_birth[1], 0.10 * w);
    real_t new_fitness2 = logcauchyv(model.m_birth[2], 0.10 * w);
    // force the fitnesses to be positive (negative values are meaningless in this context)
    new_fitness1 *= 1 - 2 * (new_fitness1 < 0);
    new_fitness2 *= 1 - 2 * (new_fitness2 < 0);
    // inipop is unidentifiable:
    real_t new_inipop = model.m_initial_pops[0];

    return instantiate_model(new_rloh, new_mu, new_fitness1,
                             new_fitness2, new_inipop);
}

/* MINIMISATION METHODS (should be in their own compile unit/library)
 * Valid method_min methods have the signature:
 *  Model method_min(std::function<real_t(Model &model)> objective,
 *                   Model  initial_guess);
 * The first argument is a closure to this object -- this allows the other
 * parameters of the log-likelihood objective function to be "bound" in the
 * context. We then vary with regards to the (Model) guess to find the
 * best_guess.
 *
 * The signature of the minimisation method is independent of the specific
 * algorithm used.
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

    // Initialise variables:
    Model best_guess = initial_guess;
    double best_score = objective(best_guess);

    // Width of neighbourhood:
    double w = log(2); // initial width for log-cauchy variates

    printf("Initial likelihood:\n-log L = %g\n", best_score);

    // Simulated annealing process:
    double Temp = best_score;      // Initial temperature
    unsigned int iter = 0;  // count iterations
    unsigned int reheating_tries = 12;

    while (Temp > Tmin) {
        Model new_guess = get_neighbour(best_guess, w);

        double score_new = objective(new_guess);

        if (std::isnan(score_new)) continue;

        double delta_score = score_new - best_score;

        if ((delta_score < 0) || ((rand() / (RAND_MAX + 1.0)) < exp(-delta_score / Temp))) {
            // Update best guess:
            best_guess = new_guess;
            best_score = score_new;

            /* Shrink the width: */
            w *= 1.0 - smoothing_factor;
            w += smoothing_factor * min_width;
        }

        Temp *= delta;
        ++iter;

        if (Temp <= Tmin) {
            if (reheating_tries == 0) {
                break;
            } else {
                /* reheat: */
                Temp = 1.0;
                w = log(2);
                --reheating_tries;
            }
        }
    }

    printf("System fully cooled after %d iterations\n", iter);
    printf("-log L = %.8g\n", best_score);
    return best_guess;
}

Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess, int resolution) {
    // the brute force method samples every value in a cuboid neighbourhood of
    // the initial guess (in log space) and returns the lowest candidate.
    const double length = 10.0f; // side length of the cube: +-sqrt(this)
    // resolution is the number of taps along the side of the cube

    // Initialise parameters
    Model best_guess = initial_guess;
    std::vector<real_t> params = model_params_pure(best_guess);
    int dim = params.size();

    // find the starting corner of the cube
    double loglength = log(length);
    std::vector<real_t> Delta(dim,0);
    for (size_t index = 0; index < dim; ++index) {
        Delta[index] = params[index] * (exp(-loglength/2) - 1);
    }

    Model start_corner = shifted_model(best_guess, Delta);
    double best_score = objective(start_corner);

    // distance between samples in log space
    real_t epsilon = loglength / resolution;
    printf("tol: %.8f\n", epsilon / 2);

    // find the volume of the giant cube
    size_t volume = 1;
    for (size_t index = 0; index < dim; ++index) {
        volume *= resolution;
    }
    // volume is now about a million if resolution is 32
    // NB: volume will be limited by the max value of size_t

    // for samples from 0 to ... volume,
    for (size_t sample = 0; sample < volume; ++sample) {
        // map this sample to a unique point in the cube
        size_t axis_tap = sample % resolution;
        size_t next_tap = sample / resolution;
        // we now have the tap (along the edge of the cube) for the first axis.
        std::vector<real_t> sample_coord(dim,0);
        for (size_t axis = 0; axis < dim; ++axis) {
            sample_coord[axis] = params[axis] * (exp(axis_tap * epsilon) - 1);
            // Repeatedly divide by resolution to get the next one:
            axis_tap = next_tap % resolution;
            next_tap = next_tap / resolution;
        }
        // sample_coord now contains the coordinates of the new guess relative
        // to the starting corner (conceptually at [0,0,0,...]).

        Model new_guess = shifted_model(start_corner, sample_coord);
        double new_score = objective(new_guess);
        if (new_score < best_score) {
            best_guess = new_guess;
            best_score = new_score;
        }

        // and just keep going until the end. don't break, check them all.
    }

    printf("Sampling complete after %d samples\n", volume);
    printf("-log L = %.8g\n", objective(best_guess));

    return best_guess;
}

Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess) {
    return brute_force_min(objective, initial_guess, 4);
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

Eigen::MatrixXd gradient_log(std::function<real_t(Model&)> objective,
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
    std::vector<real_t> Theta = model_params_pure(point);

    for (int axis = 0; axis < dim; ++axis) {
        /* Compute the derivative along each axis using a finite difference
           scheme and the weights given above: */
        double diff = 0;
        for (int tap = 0, step = -radius; tap < weights.size(); ++tap, ++step) {
            std::vector<real_t> Delta(dim, 0);
            Delta[axis] = epsilon * Theta[axis] * step;
            Model dmodel = shifted_model(point, Delta);
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
    std::vector<real_t> Theta = model_params_pure(point);

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
                Model dmodel = shifted_model(point, Delta);
                diff += tap[0] * objective(dmodel) / epsilon;
            }
            // Scale diff appropriately to get second derivative:
            diff /= epsilon;

            Hessian(x,y) = diff;
        }
    }

    //return Hessian; //?
    // TODO reconsider what units to output the hessian in
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
    double best_score = objective(best_guess);

    int dim = 4;
    double learning_rate = 1e-4;
    double tolerance = 1e-3; // minimum "speed": % change tolerated in parameters

    // while the gradient is too large:
    while (1) {
        // compute the gradient at the current best guess
        Eigen::MatrixXd gradient = gradient_log(objective, best_guess);

        std::vector<real_t> Theta = model_params_pure(best_guess);

        // update the current best guess by -learning_rate * gradient
        std::vector<real_t> Delta(dim, 0);
        for (int axis = 0; axis < dim; ++axis) {
            Delta[axis] = Theta[axis] * (exp(-learning_rate * gradient(axis)) - 1);
        }

        Model new_guess = shifted_model(best_guess, Delta);

        // compute the change in score
        double dscore = objective(new_guess) - best_score;
        // the score should always decrease: if it increases we have overshot
        // the minimum, and should try again with a smaller learning rate

        if (dscore > 0) {
            learning_rate *= 0.5;
            continue;
        }

        best_guess = new_guess;
        best_score = objective(best_guess);

        // otherwise, we want the percentage-change in all coefficients to be
        // lower than the tolerance:

        Eigen::MatrixXd Jacobian(dim,dim);
        double change = 0;
        for (int i = 0; i < dim; ++i) {
            Jacobian(i,i) = 1.0 / Theta[i];
            double dx = Delta[i] / Theta[i];
            change += dx * dx;
        }

        change = std::sqrt(change);

        // if the gradient is small enough, exit:
        if (change < tolerance) break;
        if (std::isnan(change)) break;
    }

    printf("-log L = %.8g\n", objective(best_guess));

    return best_guess;
}

Model mixed_min(std::function<real_t(Model& model)> objective,
                Model initial_guess) {
    // Minimisation method that follows one pass of brute-force search with one
    // pass of gradient descent
    Model better_guess = brute_force_min(objective, initial_guess);
    return gradient_min(objective, better_guess);
}

Model mixed_min_8(std::function<real_t(Model& model)> objective,
                Model initial_guess) {
    Model better_guess = brute_force_min(objective, initial_guess, 8);
    return gradient_min(objective, better_guess);
}

Model mixed_min_16(std::function<real_t(Model& model)> objective,
                Model initial_guess) {
    Model better_guess = brute_force_min(objective, initial_guess, 16);
    return gradient_min(objective, better_guess);
}

Model skip_minimisation(std::function<real_t(Model& model)> objective,
                        Model initial_guess) {
    // don't bother minimising, just return the initial guess.
    return initial_guess;
}

// end of minimisation methods

// Methods to serialise/output data:

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
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence,
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
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence) {
    // save with a default filename:
    save_histogram(binwidth, max_age, reference_pop, end_nodes, incidence,
                   "syntheticdata_hist.csv");
}

std::vector<size_t> parseStringToVector(const std::string& str) {
    std::vector<size_t> result;
    std::stringstream ss(str);
    char c;

    while (ss.good()) {
        std::string numStr;
        ss >> c;

        if (c == '[') continue; // Skip '['
        if (c == ']') break; // Stop at ']'

        while (c != ',' && c != ']') {
            numStr += c;
            ss >> c;
        }

        result.push_back(std::stoul(numStr));
    }

    return result;
}

void load_histogram(real_t& binwidth, real_t& max_age, size_t& reference_pop,
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence,
                    std::string filename) {
    std::ifstream histogram;
    histogram.open(filename);

    auto findset = [](std::string line, const char* match, auto& varbl) {
        int posn = line.find(match);
        if (posn != std::string::npos) {
            varbl = std::stod(line.substr(line.find("=") + 1));
            // TODO check stoi/stod behaviour is correct
        }
    };

    std::string line;
    // read the first three lines of the file and use them to set the three
    // scalar variables:
    for (size_t linenr = 0; linenr < 3; ++linenr) {
        std::getline(histogram, line);
        findset(line, "bin width", binwidth);
        findset(line, "max age", max_age);
        findset(line, "ref population", reference_pop);
    }
    std::cout << "bin width = " << binwidth << std::endl;
    std::cout << "max age = " << max_age << std::endl;
    std::cout << "ref population = " << reference_pop << std::endl;

    // then read in bars from histogram:
    for (auto &end_node : end_nodes) {
        std::getline(histogram, line);
        std::cout << end_node << ": ";
        std::cout << line << std::endl;
        incidence[end_node] = parseStringToVector(line);
    }

    histogram.close();
}

Histogram_t jackknife_incidence(size_t index, const Histogram_t& histogram,
                                std::vector<size_t> end_nodes) {
    // delete the index-th entry in the histogram
    Histogram_t new_incidence = histogram;
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
    // if we get here something has gone wrong and there are not enough entries
    // in the histogram.
    return new_incidence;
}

void resample_incidence(
    Model (*method_min)(std::function<real_t(Model&)>, Model),
    const Histogram_t *incidence,
    size_t reference_pop, size_t start, size_t end, std::vector<size_t> end_nodes,
    real_t binwidth, Model *initial_guess, std::vector<Model> *resampled_estimates) {

    for (unsigned tries = start; tries < end; ++tries) {
        // Resample the incidence:
        Histogram_t resampled_incidence;
        resampled_incidence = jackknife_incidence(tries, *incidence, end_nodes);

        std::cout << "Resampling..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, resampled_incidence);
            // the histogram version
        };

        // use a chosen method_min method to minimise the objective:
        Model best_guess = method_min(objective, *initial_guess);
        // Minimisation now complete. Print guess:
        std::cout << "Best guesses:" << std::endl;
        print_model(best_guess);
        // Annealing now complete. Store guessed model parameters:
        resampled_estimates->push_back(best_guess);
    }
}

void jackknife_and_save(Model (*method_min)(std::function<real_t(Model&)>, Model),
                        Histogram_t &incidence,
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
                                            method_min,
                                            &incidence, reference_pop,
                                            start, end, end_nodes, binwidth,
                                            &initial_guess, &results[thread]));
    }

    // run runs_per_thr + remainder in this, the parent thread:
    resample_incidence(method_min, &incidence, reference_pop,
                       reference_pop - runs_per_thr - remainder,
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

real_t save_data_compute_maximum(Epidata_t &all_times) {
    // compute the maximum age and save the all_times data:
    real_t max_age = 0;

    std::ofstream fakedata;
    fakedata.open("syntheticdata_raw.csv");
    for (auto& entry : all_times) {
        fakedata << entry.first << "," << entry.second << "," << std::endl;
        max_age += (entry.first > max_age) * (entry.first - max_age);
    }
    fakedata.close();

    return max_age;
}

void draw_level_sets(std::function<real_t(Model &model)> objective,
                     Model point, int x_axis, int y_axis,
                     real_t x_range=10.0f, real_t y_range=10.0f,
                     size_t lines=16) {
    using Eigen::all;
    // Use this function to draw realistic non-ellipsoidal likelihood intervals.
    // We choose two axes, and want to know about the projection of a surface
    // objective=const onto this plane. Along this surface, because we
    // solve objective = const., d obj / d marginal_variables = 0.
    std::vector<real_t> params = model_params_pure(point);

    // Create the paper:
    std::ofstream drawing;
    char filename[25];
    std::snprintf(filename, sizeof(filename), "raster_intervals_%d_%d.csv", x_axis, y_axis);
    drawing.open(filename);

    drawing << "x,y,z*," << std::endl;

    // Copy the model object:
    Model sample_point = point;

    // pick a pixel, and instead of sampling the objective function here, find
    // the minimum in the normal directions.
    //    we try to minimise H along the normal directions
    //    because when H is minimised dH/dn = 0, and this is one of the
    //    conditions of the projection.
    // Store the normal directions for later (the ones that are not chosen):
    std::vector<int> normal_dirs{0,1,2,3};
    for (auto dir = normal_dirs.begin(); dir != normal_dirs.end();) {
        if ((*dir == x_axis) || (*dir == y_axis)) {
            normal_dirs.erase(dir);
            continue;
        } else {
            ++dir;
        }
    }

    for (int x_tap = 0; x_tap < lines; ++x_tap) {
        for (int y_tap = 0; y_tap < lines; ++y_tap) {
            // compute new parameters in log grid:
            real_t x_value = params[x_axis] / x_range;
            real_t y_value = params[y_axis] / y_range;
            x_value *= exp(log(x_range) * 2 * (real_t)x_tap / (real_t)lines);
            y_value *= exp(log(y_range) * 2 * (real_t)y_tap / (real_t)lines);

            // store current parameters:
            std::vector<real_t> sample_params = model_params_pure(point);

            // compute differences:
            real_t dx, dy;
            dx = x_value - sample_params[x_axis];
            dy = y_value - sample_params[y_axis];
            std::vector<real_t> Delta(4,0);
            Delta[x_axis] = dx, Delta[y_axis] = dy;

            // Move the sample point to have the right x- and y-values:
            sample_point = shifted_model(point, Delta);

            real_t z_value = objective(sample_point);

            // instead of just sampling the objective at this point, minimise it
            // with respect to the normal directions (the parameters we are
            // ignoring)
            int dim = 4;
            double learning_rate = 1;
            // Cap the number of newton's method iterations at 4
            int max_iter = 4;

            for (int iter = 0; iter < max_iter; ++iter) {
                // Minimise it with Newton's method
                Eigen::MatrixXd gradient = gradient_log(objective, sample_point);
                Eigen::MatrixXd normal_gradient = gradient(normal_dirs, all);

                // Get Hessian for newtonian steps:
                Eigen::MatrixXd hessian = compute_hessian(objective, sample_point);
                Eigen::MatrixXd normal_hessian = hessian(normal_dirs, normal_dirs);

                // Compute Newtonian step:
                Eigen::MatrixXd newt_step = normal_hessian.inverse() * normal_gradient;
                newt_step = newt_step * -1 * learning_rate;

                // Update the guess:
                std::vector<real_t> Step(dim, 0);
                sample_params = model_params_pure(sample_point);
                for (int axis = 0; axis < 2; ++axis) {
                    Step[normal_dirs[axis]] = newt_step(axis);
                    Step[normal_dirs[axis]] /= sample_params[normal_dirs[axis]];
                }

                Model new_point = shifted_model(sample_point, Step);

                // quit when we start to go uphill:
                real_t new_z_value = objective(new_point);

                if (new_z_value > z_value) break;
                if (std::isnan(new_z_value)) break;

                sample_point = new_point;
                z_value = new_z_value;
            }

            // write out values:
            drawing << x_value << "," << y_value << "," << z_value << ",";
            drawing << std::endl;
        }
    }

    drawing.close();
}

void draw_3dsurface(std::function<real_t(Model &model)> objective,
                    Model origin, int x_axis, int y_axis,
                    real_t x_range=10.0f, real_t y_range=10.0f,
                    size_t lines=16) {
    /* This function samples the objective function in a neighbourhood around
     * the origin point. It then outputs a CSV of the form
     * x, y, z,
     * where z = objective(x, y)
     */
    std::ofstream drawing;
    char filename[14];
    std::snprintf(filename, sizeof(filename), "mesh_%d_%d.csv", x_axis, y_axis);
    drawing.open(filename);
    drawing << "x,y,z," << std::endl;
    // We will need to translate between different coordinate axes and Model
    // parameters, so here is a vector of parameters:
    std::vector<real_t> orig_params = model_params_pure(origin);
    /* ...
    * sample points on a logarithmic scale from origin[x]/x_range to
    * origin[x]*x_range and the same for y.
    */
    for (int x_tap = 0; x_tap < lines; ++x_tap) {
        for (int y_tap = 0; y_tap < lines; ++y_tap) {
            // compute new parameters in log grid:
            real_t x_value = orig_params[x_axis] / x_range;
            real_t y_value = orig_params[y_axis] / y_range;
            x_value *= exp(log(x_range) * 2 * (real_t)x_tap / (real_t)lines);
            y_value *= exp(log(y_range) * 2 * (real_t)y_tap / (real_t)lines);

            // compute differences:
            real_t dx, dy;
            dx = x_value - orig_params[x_axis];
            dy = y_value - orig_params[y_axis];
            std::vector<real_t> Delta(4,0);
            Delta[x_axis] = dx, Delta[y_axis] = dy;

            Model sample_point = shifted_model(origin, Delta);

            real_t z_value = objective(sample_point);

            // write out values:
            drawing << x_value << "," << y_value << "," << z_value << ",";
            drawing << std::endl;
        }
    }
    drawing.close();
}

// The bounding box has a centre, side dimensions, and resolutions.
// The axes are mu, rloh, s1, and the colour will be given by s1 and the
// likelihood.
struct BoundingBox {
    unsigned int resolution[3];
    double centre[3];
    double dims[3];
    double s2;
    size_t size;
};

// TODO: reference to float vector safer than raw pointer?
void sample_voxel_cube(float* buffer,
                       std::function<real_t(Model &model)> objective,
                       struct BoundingBox bounding_box,
                       real_t offset, real_t s2max,
                       size_t first_sample, size_t last_sample) {
    // gamma compression:
    float compression = 10.0 / bounding_box.size; // FIXME ad hoc

    for (size_t sample = first_sample; sample < last_sample; ++sample) {
        // compute row, col, and lyr
        size_t next;
        size_t row = sample % bounding_box.resolution[0];
        next = sample / bounding_box.resolution[0];
        size_t col = next % bounding_box.resolution[1];
        next = next / bounding_box.resolution[1];
        size_t lyr = next % bounding_box.resolution[2];
        next = next / bounding_box.resolution[2];

        // set position in box
        double offset_mu = row * 1.0 / bounding_box.resolution[0];
        offset_mu -= 0.5;
        real_t mu = bounding_box.centre[0]; 
        mu *= exp(log(bounding_box.dims[0]) * offset_mu);

        double offset_rloh = col * 1.0 / bounding_box.resolution[1];
        offset_rloh -= 0.5;
        real_t rloh = bounding_box.centre[1]; 
        rloh *= exp(log(bounding_box.dims[1]) * offset_rloh);

        double offset_s1 = lyr * 1.0 / bounding_box.resolution[2];
        offset_s1 -= 0.5;
        offset_s1 *= 2.0 * bounding_box.dims[2];
        real_t s1 = bounding_box.centre[2];
        s1 *= (1.0 + offset_s1);

        // Marginalise/integrate over s2:
        real_t slices = 16; // TODO I don't see a big difference between 8 and 64
        real_t ds2 = s2max / slices;
        float voxel[3] = {0, 0, 0};
        for (real_t s2 = 0.0; s2 < s2max; s2 += ds2) {
            // Associate a floating point colour with this s2 value:
            float theta = (4 * M_PI / 3) * s2 / s2max;
            // swatch 1:
            /*float colour[3] = {0.66f -0.33f * cosf(theta),
                               0.66f +0.17f * cosf(theta) -0.29f * sinf(theta),
                               0.66f +0.17f * cosf(theta) +0.29f * sinf(theta)};
                               */
            // swatch 2:
            float colour[3] = {0.5f -0.408f * cosf(theta),
                               0.5f +0.204f * cosf(theta) -0.353f * sinf(theta),
                               0.5f +0.204f * cosf(theta) +0.353f * sinf(theta)};
            // swatch 3:
            /*float colour[3] = {0.5f -0.392f * cosf(theta) +0.276f * sinf(theta),
                               0.5f +0.087f * cosf(theta) -0.122f * sinf(theta), 
                               0.5f +0.298f * cosf(theta) +0.398f * sinf(theta)};
                               */

            // get a normalised value for likelihood, 
            // which is exp(-log(L)).
            Model sample_point = instantiate_model(rloh, mu, s1, s2, 1e6);
            float dloglike = offset - objective(sample_point);
            // Apply gamma compression to squash the value into a visible range:
            dloglike *= compression;
            float y_value = exp(dloglike);

            // Nice the output: set NaN values to zero
            // Values are also invalid iff they would cause an overflow in exp:
            // 127 * log(2) = 88 (exponent should not go out of range)

            if (!std::isnormal(y_value) || dloglike > 88) {
                y_value = 0.0f;
            }
            // increment integral and convert value to floating point colour:
            for (int ch = 0; ch < 3; ++ch) {
                voxel[ch] += y_value * colour[ch] * ds2 / s2max;
                // = \frac{1}{s2max} \int_{s2=0}^s2max ... ds2
            }
        }

        // copy colour to buffer:
        std::memcpy(buffer + 3 * sample, voxel, sizeof(voxel));
    }
}

void render_voxel_cube(std::function<real_t(Model &model)> objective,
                       struct BoundingBox bounding_box,
                       std::string voxel_file,
                       size_t num_child_threads) {
    // Serialise the objective function to a voxel cube file
    real_t s2max = 0.20; // WARNING values too high result in garbage
    Model centre_of_box = instantiate_model(bounding_box.centre[1],
                                            bounding_box.centre[0],
                                            bounding_box.centre[2],
                                            bounding_box.s2,
                                            1e6);
    // likelihood functions are only unique up to a scalar constant. We use this
    // to wrap the likelihood function into a representable range by offsetting
    // the log-likelihood by a constant before converting it to a probability
    // with exp(- -loglikelihood): i.e. exp(offset - -loglikelihood)
    // note the value of the objective in the centre, call this the offset:
    real_t offset = objective(centre_of_box);
    // and try to normalise the colours:
    //real_t offset = +0.5 * bounding_box.size * log(bounding_box.size);
    printf("size: %u\n", bounding_box.size);
    printf("normalisation offset (= L(centre)): %g\n", offset);

    // get the total number of samples to take:
    size_t volume_samples = 1;
    for (int index = 0; index < 3; ++index) {
        volume_samples *= bounding_box.resolution[index];
    }
    // create a buffer to store results:
    float* buffer; // vectorised buffer for storing 3D stuff
    size_t buffer_size = volume_samples * 3 * sizeof(float);
    buffer = (float*)malloc(buffer_size);

    // Vectorise the sampling so that we can easily parallelise
    // multi thread: but divide up the work so 0 child threads recovers default
    // behaviour. When num_child_threads == 0, we will skip the loop where we
    // dispatch worker threads, and should run all the samples on the main
    // thread. child_samples should work out to volume_samples, and remainder
    // should work out to 0.
    // Get some number of samples from each child thread, and the remainder from
    // the parent thread:
    size_t child_samples = volume_samples / (num_child_threads + 1);
    size_t remainder     = volume_samples % (num_child_threads + 1);

    std::vector<std::thread> child_threads(0);
    for (size_t thrd = 0; thrd < num_child_threads; ++thrd) {
        size_t start, end;
        start = thrd * child_samples;
        end = start + child_samples;

        child_threads.push_back(std::thread(sample_voxel_cube,
                                            buffer, objective, bounding_box,
                                            offset, s2max, start, end));
    }

    // run the remaining threads in the parent thread:
    sample_voxel_cube(buffer, objective, bounding_box, offset, s2max, 
                      volume_samples - child_samples - remainder, 
                      volume_samples);

    // Wait for child threads to finish
    for (auto& thread : child_threads) {
        thread.join();
    }

    // write out results in buffer to file:
    std::ofstream voxelfile;
    voxelfile.open(voxel_file, std::ios::out|std::ios::binary);

    // Write header for voxel cube file:
    const char* header = "Voxel\n";
    voxelfile.write((char*)header, 6 * sizeof(char));

    voxelfile.write((char*)(bounding_box.resolution), 3 * sizeof(unsigned));

    // Write buffer:
    voxelfile.write((char*)(buffer), buffer_size);

    voxelfile.close();

    free(buffer);
    return;
}

struct Estimate {
    Model best_guess;
    Eigen::MatrixXd Hessian;
    Estimate();
    Estimate(const Model& guess) : best_guess(guess) {}
    Estimate(const Estimate& est) : best_guess(est.best_guess), Hessian(est.Hessian) {}
};

void print_best_guess(Estimate estimate) {
    std::cout << "Best guesses:" << std::endl;
    print_model(estimate.best_guess);

    std::cout << "H = " << std::endl;
    std::cout << "[rloh, mu, s1, s2]" << std::endl;
    std::cout << estimate.Hessian << std::endl;

    printf("\n");

    // Invert Hessian to get covariance matrix:
    std::cout << "cov = H^-1 = " << std::endl;
    std::cout << estimate.Hessian.inverse() << std::endl;

    // Compute confidence intervals:
    std::cout << "Standard deviations:" << std::endl;
    std::vector<real_t> Theta = model_params_pure(estimate.best_guess);
    for (int param = 0; param < Theta.size(); ++param) {
        std::cout << Theta[param] << " +/- ";
        std::cout << sqrt(estimate.Hessian.inverse()(param, param));
        std::cout << std::endl;
    }
}

Estimate get_estimate_germline(real_t binwidth, size_t reference_pop,
                               Histogram_t incidence,
                               Histogram_t incidence_germline,
                               Model (*method_min)(std::function<real_t(Model&)>, Model)) {
    std::cout << "--------------" << std::endl;
    std::cout << "Best estimate:" << std::endl;
    // Guess some initial model parameters:
    real_t rloh = 1e-7;
    real_t mu = 1e-8;
    real_t fitness1 = 0.02;
    real_t fitness2 = 0.02;
    real_t initialpop = 1e6;

    Model guess = instantiate_model(rloh, mu, fitness1, fitness2, initialpop);

    // Minimise the -log likelihood:
    std::function<real_t(Model&)> objective = [&](Model& model) {
        return loglikelihood_hist_both(model, binwidth,
                                       reference_pop, incidence,
                                       incidence_germline);
    };
    // this should now work with germline data

    std::cout << "Minimising likelihood..." << std::endl;
    Model best_guess = method_min(objective, guess);
    // Get and return Hessian:
    Eigen::MatrixXd Hessian = compute_hessian(objective, best_guess);

    // Annealing now complete. Print guess:
    Estimate estimate(best_guess);
    estimate.Hessian = Hessian;
    print_best_guess(estimate);

    return estimate;
}

void generate_histogram_germline(Model &ground_truth,
                        Model &ground_truth_germline,
                        size_t seed, size_t dataset_size,
                        const std::vector<size_t> &end_nodes,
                        real_t &max_age, size_t &reference_pop,
                        real_t &binwidth,
                        Histogram_t &incidence,
                        Histogram_t &incidence_germline) {
    Epidata_t all_times =
        generate_dataset(ground_truth, seed, dataset_size / 2);
    // to include germline mutations:
    // vary the initial populations so that instead of node 0 the cells are
    // initially on node 1
    Epidata_t all_times_germline =
        generate_dataset(ground_truth_germline, seed, dataset_size / 2);

    std::cout << "Done. Saving..." << std::endl;

    // compute maximum age and save the all_times data:
    max_age = save_data_compute_maximum(all_times);

    // Convert age data to histogram:
    reference_pop = all_times.size(); /* NB: with germline this is 1/2
    the previous value */
    binwidth = max_age / (2 * pow(reference_pop, 0.4)); // years
    std::cout << "max age = " << max_age;
    std::cout << "\n bin width = " << binwidth << std::endl;

    for (auto &end_node : end_nodes) {
        incidence[end_node] = convert_to_histogram(all_times, binwidth,
                              end_node);
        incidence_germline[end_node] = convert_to_histogram(
                                           all_times_germline, binwidth,
                                           end_node);
    }

    // Save the histogram:
    save_histogram(binwidth, max_age, reference_pop, end_nodes,
                   incidence, "syntheticdata_hist_sporadic.csv");
    save_histogram(binwidth, max_age, reference_pop, end_nodes,
                   incidence_germline, "syntheticdata_hist_germline.csv");
}

void guess_parameters_germline(Model &ground_truth, GuesserConfig options,
                               Model (*method_min)(std::function<real_t(Model&)>, Model)) {
    // The main test harness for statistical inference:
    // * Generate data from a clinical study
    // * Minimise the -log likelihood function
    // * Resample the estimates (if appropriate)
    // this version includes germline cases, the clinical study has a 50:50 mix of
    // sporadic cases and cases with germline alterations
    Model ground_truth_germline = ground_truth;
    // For germline cases, move the initial cells to the state with one
    // mutation:
    ground_truth_germline.m_initial_pops[1] = ground_truth.m_initial_pops[0];
    ground_truth_germline.m_initial_pops[0] = 0;
    // fitnesses 1 and 2 = 0.05, 0.03, // these are set to zero: discuss with
    // Miriam and David
    ground_truth_germline.m_birth[1] = 0;
    ground_truth_germline.m_birth[2] = 0;

    size_t seed = options.seed;
    size_t dataset_size = options.dataset_size;

    printf("Ground truth (germline):\n");
    print_model(ground_truth_germline);

    std::vector<size_t> end_nodes = {3,4};
    Histogram_t incidence_sporadic, incidence_germline;
    real_t binwidth, max_age;
    size_t reference_pop;

    // TODO loadable histograms for germline studies
    if (!options.histogram_file.empty() &&
        !options.germline_histogram_file.empty()) {
        load_histogram(binwidth, max_age, reference_pop, end_nodes,
                       incidence_sporadic, options.histogram_file);
        load_histogram(binwidth, max_age, reference_pop, end_nodes,
                       incidence_germline, options.germline_histogram_file);
    } else {
        generate_histogram_germline(ground_truth, ground_truth_germline,
                           seed, dataset_size, end_nodes, max_age, reference_pop,
                           binwidth, incidence_sporadic, incidence_germline);
    }

    // Get the main un-resampled estimate with simulated annealing:
    printf("Target likelihood:\n-log L = %g\n",
           loglikelihood_hist_both(ground_truth, binwidth, reference_pop,
                                   incidence_sporadic, incidence_germline));

    Estimate estimate = get_estimate_germline(binwidth, reference_pop,
                        incidence_sporadic, incidence_germline, method_min);

    // Annealing now complete. Print guess:
    print_best_guess(estimate);

    // Save the model to a file:
    if (!options.estimate_file.empty()) {
        clonal_expansion::save_model(estimate.best_guess, options.estimate_file);
    }

    if (options.resample_after) {
        std::cout << "resampling not yet supported for germline studies";
        std::cout << std::endl;
    }

    // LIKELIHOOD FUNCTION VISUALISATION (optional)

    // Draw level sets:
    if (options.level_sets) {
        int dim = 4;
        std::cout << "Tracing level sets..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, incidence_sporadic,
                                           incidence_germline);
        };
        // choose a starting point from the best guess and local hessian:
        Model start_point = estimate.best_guess;

        // for each combination of parameters:
        for (int x_axis = 0; x_axis < dim; ++x_axis) {
            for (int y_axis = x_axis + 1; y_axis < dim; ++y_axis) {
                draw_level_sets(objective, start_point, x_axis, y_axis);
            }
        }
    }

    // Draw 3d voxel cube for likelihood function:
    if (!options.voxel_file.empty()) {
        struct BoundingBox bbox = {
            // For now, just assume the same resolution along each axis:
            .resolution = {options.voxel_res, options.voxel_res,
                          options.voxel_res},
            .centre = {estimate.best_guess.m_migr[0][1],
                       estimate.best_guess.m_migr[0][2],
                       estimate.best_guess.m_birth[1]
                     },
            // TODO pass these in on command line
            // Also, use logarithmic scales for mu and rloh and linear scales for
            // s1 and s2.
            .dims = {options.mesh_x_range, options.mesh_y_range, 1},
            .s2 = ground_truth.m_birth[2],
            .size = reference_pop
        };

        std::cout << "Rendering 3D voxel data..." << std::endl;

        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth, reference_pop, 
                                           incidence_sporadic, incidence_germline);
        };

        render_voxel_cube(objective, bbox, options.voxel_file, 
                          options.num_child_threads);
    }
}

Estimate get_estimate(real_t binwidth, size_t reference_pop,
                      Histogram_t incidence,
                      Model (*method_min)(std::function<real_t(Model&)>, Model)) {
    std::cout << "--------------" << std::endl;
    std::cout << "Best estimate:" << std::endl;
    // Guess some initial model parameters:
    real_t rloh = 1e-6;
    real_t mu = 1e-7;
    real_t fitness1 = 0.03;
    real_t fitness2 = 0.03;
    real_t initialpop = 1e6;

    Model guess = instantiate_model(rloh, mu, fitness1, fitness2, initialpop);

    // Minimise the -log likelihood:
    std::function<real_t(Model&)> objective = [&](Model& model) {
        return loglikelihood_hist_both(model, binwidth,
                                       reference_pop, incidence);
    };

    std::cout << "Starting minimisation..." << std::endl;
    Model best_guess = method_min(objective, guess);
    // Get and return Hessian:
    Eigen::MatrixXd Hessian = compute_hessian(objective, best_guess);

    Estimate estimate(best_guess);
    estimate.Hessian = Hessian;

    return estimate;
}

void generate_histogram(Model &ground_truth, size_t seed, size_t dataset_size,
                        const std::vector<size_t> &end_nodes,
                        real_t &max_age, size_t &reference_pop,
                        real_t &binwidth, Histogram_t &incidence) {
    // simulate some data that could be collected by a longitudinal study:
    std::cout << "Generating synthetic dataset..." << std::endl;
    // sporadic cases only (default)
    Epidata_t all_times = generate_dataset(ground_truth, seed, dataset_size);

    std::cout << "Done. Saving..." << std::endl;

    // compute maximum age and save the all_times data:
    max_age = save_data_compute_maximum(all_times);

    // Convert age data to histogram:
    reference_pop = all_times.size();
    binwidth = max_age / (2 * pow(reference_pop, 0.4)); // years
    //real_t binwidth = 10.0;
    std::cout << "max age = " << max_age;
    std::cout << "\n bin width = " << binwidth << std::endl;

    // set the incidence for each karyotype:
    for (auto &end_node : end_nodes) {
        incidence[end_node] = convert_to_histogram(all_times, binwidth,
                              end_node);
    }

    // Save the histogram:
    save_histogram(binwidth, max_age, reference_pop, end_nodes, incidence);
}

void guess_parameters(Model &ground_truth, GuesserConfig options,
                      Model (*method_min)(std::function<real_t(Model&)>, Model)) {
    // The main test harness for statistical inference:
    // * Generate data from a clinical study
    // * Minimise the -log likelihood function
    // * Resample the estimates (if selected)
    printf("Ground truth:\n");
    print_model(ground_truth);
    size_t seed = options.seed;
    size_t dataset_size = options.dataset_size;

    std::vector<size_t> end_nodes = {3,4};
    Histogram_t incidence;
    real_t binwidth, max_age;
    size_t reference_pop;

    // loadable histograms:
    if (!options.histogram_file.empty()) {
        std::cout << options.histogram_file << std::endl;
        load_histogram(binwidth, max_age, reference_pop, end_nodes, incidence,
                       options.histogram_file);
    } else {
        generate_histogram(ground_truth, seed, dataset_size, end_nodes,
                           max_age, reference_pop, binwidth, incidence);
    }

    printf("Target likelihood:\n-log L = %g\n",
           loglikelihood_hist_both(ground_truth, binwidth, reference_pop,
                                   incidence));

    // Get the main un-resampled best parameter estimate:
    Estimate estimate = get_estimate(binwidth, reference_pop, incidence,
                                     method_min);

    // Annealing now complete. Print guess:
    print_best_guess(estimate);

    // Save the model to a file:
    if (!options.estimate_file.empty()) {
        clonal_expansion::save_model(estimate.best_guess, options.estimate_file);
    }

    // Jack-knife resampling of incidence:
    if (options.resample_after) {
        jackknife_and_save(method_min, incidence, reference_pop, binwidth, end_nodes,
                           estimate.best_guess, options.num_child_threads);
    }

    // LIKELIHOOD FUNCTION VISUALISATION (optional)

    // Draw level sets:
    if (options.level_sets) {
        int dim = 4;
        std::cout << "Tracing level sets..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, incidence);
            // Sanity check test function:
            //real_t dx = log(model.m_migr[0][2] / ground_truth.m_migr[0][2]);
            //real_t dy = log(model.m_migr[0][1] / ground_truth.m_migr[0][1]);
            //return dx * dx + dy * dy;
        };
        // choose a starting point from the best guess and local hessian:
        Model start_point = estimate.best_guess;

        // for each combination of parameters:
        for (int x_axis = 0; x_axis < dim; ++x_axis) {
            for (int y_axis = x_axis + 1; y_axis < dim; ++y_axis) {
                draw_level_sets(objective, start_point, x_axis, y_axis,
                                options.mesh_x_range,
                                options.mesh_y_range,
                                options.mesh_lines);
            }
        }
    }

    // Draw 3d plot:
    if (options.draw_mesh) {
        int dim = 4;
        std::cout << "Sampling likelihood function..." << std::endl;
        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, incidence);
        };

        // for each combination of parameters:
        for (int x_axis = 0; x_axis < dim; ++x_axis) {
            for (int y_axis = x_axis + 1; y_axis < dim; ++y_axis) {
                draw_3dsurface(objective, estimate.best_guess, x_axis, y_axis,
                               options.mesh_x_range,
                               options.mesh_y_range,
                               options.mesh_lines);
            }
        }
    }

    // Draw 3d voxel cube:
    if (!options.voxel_file.empty()) {
        struct BoundingBox bbox = {
            // For now, just assume the same resolution along each axis:
            .resolution = {options.voxel_res, options.voxel_res,
                          options.voxel_res},
            // use ground_truth as the centre of the box:
            .centre = {ground_truth.m_migr[0][1],
                       ground_truth.m_migr[0][2],
                       ground_truth.m_birth[1]
                     },
            // TODO pass ground truth/initial guess in on command line?
            // Also, use logarithmic scales for mu and rloh and linear scales for
            // s1 and s2. DONE
            .dims = {options.mesh_x_range, options.mesh_y_range, 1},
            .s2 = ground_truth.m_birth[2],
            .size = reference_pop
        };

        std::cout << "Rendering 3D voxel data..." << std::endl;

        std::function<real_t(Model&)> objective = [&](Model& model) {
            return loglikelihood_hist_both(model, binwidth,
                                           reference_pop, incidence);
        };

        render_voxel_cube(objective, bbox, options.voxel_file, 
                          options.num_child_threads);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Call this program with:\n";
        std::cout << "./guesser --seed seed --sample_size dataset_size ...\n";
        std::cout << std::endl;
        return 1;
    }
    // read in config options:
    GuesserConfig options(argc, argv);

    srand(options.seed);

    Model ground_truth = instantiate_model(5.0e-7,
                                           5.0e-8,
                                           0.05,
                                           0.03,
                                           1e6);
    // TODO load initial guess from a file:

    // the default minimisation method is annealing:
    typedef Model (*Minimiser_t)(std::function<real_t(Model&)>, Model);
    Minimiser_t method_min = annealing_min;
    // but:
    if (options.minimise_with_gradient) method_min = gradient_min;
    if (options.minimise_brute_force)   method_min = brute_force_min;
    if (options.minimise_with_mixed)    method_min = mixed_min;
    if (options.minimise_with_mixed_8)  method_min = mixed_min_8;
    if (options.minimise_with_mixed_16) method_min = mixed_min_16;
    if (options.skip_minimisation) method_min = skip_minimisation;

    // The inference harness itself:
    void (*guessing_harness)(Model &ground_truth, GuesserConfig options,
                             Minimiser_t method_min);
    guessing_harness = guess_parameters; // default value

    if (options.include_germline)
        guessing_harness = guess_parameters_germline;

    // run the inference harness:
    guessing_harness(ground_truth, options, method_min);

    return 0;
}
