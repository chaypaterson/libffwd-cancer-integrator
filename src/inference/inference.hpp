#ifndef INFERENCE
#define INFERENCE

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

using clonal_expansion::real_t;
using clonal_expansion::Model;
using clonal_expansion::GuesserConfig;

typedef std::vector<std::pair<double, int>> Epidata_t;
typedef std::map<size_t, std::vector<size_t>> Histogram_t;

// Structures
struct Estimate {
    Model best_guess;
    Eigen::MatrixXd Hessian;
    Estimate(Model guess) : best_guess(guess) {}
};

struct BoundingBox {
    unsigned int resolution[3];
    real_t centre[3];
    real_t dims[3];
    real_t s2;
    size_t size;
};

// Statistical functions:

// return the survival probability for this node:
real_t logsurvival(Model& params, int node);

real_t loglikelihood_hist_node(Model& params, size_t node, real_t binwidth,
                               size_t ref_population,
                               const std::vector<size_t>& freqs);

// Histogram version of the -log likelihood
// Recieves a histogram of cancers with a known type (node)
// Returns a -log binomial likelihood
real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               Histogram_t histos);
// Histogram version of the -log likelihood
// Applied to a set of histograms of cancer incidence with given types
real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               Histogram_t histos_sporadic,
                               Histogram_t histos_germline);

std::vector<size_t> convert_to_histogram(const Epidata_t& all_times,
        real_t binwidth,
        size_t node);


// Gillespie algorithm functions:

Epidata_t generate_dataset(Model& model, int seed, int runs);

// Spawn a model of tumour suppressor loss
Model instantiate_model(real_t rloh, real_t mu, real_t fitness1,
                        real_t fitness2, real_t initialpop);

// return a vector of values that summarise the model parameters
std::vector<real_t> model_params_pure(Model point);

// return a vector of raw pointers to the model parameters
std::vector<volatile real_t*> model_params_raw(Model point);

// Spawn a model of tumour suppressor loss
Model instantiate_model_germline(real_t rloh, real_t mu, real_t fitness1,
                                 real_t fitness2, real_t initialpop);

Model shifted_model(Model& params, real_t drloh, real_t dmu, real_t dfitness1,
                    real_t dfitness2, real_t dinitialpop);

Model shifted_model(Model& params, std::vector<real_t> Delta);

double logcauchyv(double mode, double width);

double uniform(double mean, double width);

Model get_neighbour(Model& model, double w);

// Function reads a set of datum points and returns a set of model
Model annealing_min(std::function<real_t(Model &model)> objective,
                    Model initial_guess);
// the brute force method samples every value in a cuboid neighbourhood of
// the initial guess (in log space) and returns the lowest candidate.
Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess, int resolution);

Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess);

// Numerical analysis methods:

// Stencil for numerical differentiation:
std::vector<std::vector<double>> StencilLP16();

// Compute the gradient of the objective function at a point with respect to
// the logs of the model parameters:
Eigen::MatrixXd gradient_log(std::function<real_t(Model&)> objective,
                             Model point);
// Compute the hessian of the objective function at a point:
Eigen::MatrixXd compute_hessian(std::function<real_t(Model&)> objective,
                                Model point);
/* Minimise objective using gradient descent (NB: we are using numerical
       differentiation, not autodiff/backpropagation)
*/
Model gradient_min(std::function<real_t(Model& model)> objective,
                   Model initial_guess);
// Minimisation method that follows one pass of brute-force search with one
Model mixed_min(std::function<real_t(Model& model)> objective,
                Model initial_guess);

Model mixed_min_8(std::function<real_t(Model& model)> objective,
                Model initial_guess);

Model mixed_min_16(std::function<real_t(Model& model)> objective,
                Model initial_guess);
// don't bother minimising, just return the initial guess.
Model skip_minimisation(std::function<real_t(Model& model)> objective,
                        Model initial_guess);

// Methods to serialise/output data:

void print_model(Model &model);

void write_model_line(std::ofstream& file, Model &model);

void save_histogram(real_t& binwidth, real_t&  max_age, size_t& reference_pop,
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence,
                    std::string filename);

void save_histogram(real_t& binwidth, real_t&  max_age, size_t& reference_pop,
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence);

std::vector<size_t> parseStringToVector(const std::string& str);

void load_histogram(real_t& binwidth, real_t& max_age, size_t& reference_pop,
                    const std::vector<size_t>& end_nodes,
                    Histogram_t& incidence,
                    std::string filename);

Histogram_t jackknife_incidence(size_t index, const Histogram_t& histogram,
                                std::vector<size_t> end_nodes);

void resample_incidence(
    Model (*method_min)(std::function<real_t(Model&)>, Model),
    const Histogram_t *incidence,
    size_t reference_pop, size_t start, size_t end, std::vector<size_t> end_nodes,
    real_t binwidth, Model *initial_guess, std::vector<Model> *resampled_estimates);

void jackknife_and_save(Model (*method_min)(std::function<real_t(Model&)>, Model),
                        Histogram_t &incidence,
                        size_t reference_pop, real_t binwidth,
                        std::vector<size_t> end_nodes, Model &initial_guess,
                        size_t n_child_threads = 0);

// compute the maximum age and save the all_times data:
real_t save_data_compute_maximum(Epidata_t &all_times);

void draw_level_sets(std::function<real_t(Model &model)> objective,
                     Model point, int x_axis, int y_axis,
                     real_t x_range=10.0f, real_t y_range=10.0f,
                     size_t lines=16);

/* This function samples the objective function in a neighbourhood around
 * the origin point. It then outputs a CSV of the form
 * x, y, z,
 * where z = objective(x, y)
 */
void draw_3dsurface(std::function<real_t(Model &model)> objective,
                    Model origin, int x_axis, int y_axis,
                    real_t x_range=10.0f, real_t y_range=10.0f,
                    size_t lines=16);

void sample_voxel_cube(float* buffer,
                       std::function<real_t(Model &model)> objective,
                       struct BoundingBox bounding_box,
                       real_t offset, real_t s2max,
                       size_t first_sample, size_t last_sample);

// Serialise the objective function to a voxel cube file
void render_voxel_cube(std::function<real_t(Model &model)> objective,
                       struct BoundingBox bounding_box,
                       std::string voxel_file,
                       size_t num_child_threads);

void print_best_guess(Estimate estimate);

Estimate get_estimate_germline(real_t binwidth, size_t reference_pop,
                               Histogram_t incidence,
                               Histogram_t incidence_germline,
                               Model (*method_min)(std::function<real_t(Model&)>, Model));

void generate_histogram_germline(Model &ground_truth,
                        Model &ground_truth_germline,
                        size_t seed, size_t dataset_size,
                        const std::vector<size_t> &end_nodes,
                        real_t &max_age, size_t &reference_pop,
                        real_t &binwidth,
                        Histogram_t &incidence,
                        Histogram_t &incidence_germline);

// The main test harness for statistical inference:
// * Generate data from a clinical study
// * Minimise the -log likelihood function
// * Resample the estimates (if appropriate)
// this version includes germline cases, the clinical study has a 50:50 mix of
// sporadic cases and cases with germline alterations
void guess_parameters_germline(Model &ground_truth, GuesserConfig options,
                               Model (*method_min)(std::function<real_t(Model&)>, Model));

Estimate get_estimate(real_t binwidth, size_t reference_pop,
                      Histogram_t incidence,
                      Model (*method_min)(std::function<real_t(Model&)>, Model));

void generate_histogram(Model &ground_truth, size_t seed, size_t dataset_size,
                        const std::vector<size_t> &end_nodes,
                        real_t &max_age, size_t &reference_pop,
                        real_t &binwidth, Histogram_t &incidence);
// The main test harness for statistical inference:
// * Generate data from a clinical study
// * Minimise the -log likelihood function
// * Resample the estimates (if selected)
void guess_parameters(Model &ground_truth, GuesserConfig options,
                      Model (*method_min)(std::function<real_t(Model&)>, Model));

#endif
