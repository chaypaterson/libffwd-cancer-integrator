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

real_t logsurvival(Model& params, int node);

real_t loglikelihood_hist_node(Model& params, size_t node, real_t binwidth,
                               size_t ref_population,
                               const std::vector<size_t>& freqs);

real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               Histogram_t histos);

real_t loglikelihood_hist_both(Model& params, real_t binwidth,
                               size_t ref_population,
                               Histogram_t histos_sporadic,
                               Histogram_t histos_germline);

std::vector<size_t> convert_to_histogram(const Epidata_t& all_times,
        real_t binwidth,
        size_t node);

Epidata_t generate_dataset(Model& model, int seed, int runs);

Model instantiate_model(real_t rloh, real_t mu, real_t fitness1,
                        real_t fitness2, real_t initialpop);

std::vector<real_t> model_params_pure(Model point);

std::vector<volatile real_t*> model_params_raw(Model point);

Model instantiate_model_germline(real_t rloh, real_t mu, real_t fitness1,
                                 real_t fitness2, real_t initialpop);

Model shifted_model(Model& params, real_t drloh, real_t dmu, real_t dfitness1,
                    real_t dfitness2, real_t dinitialpop);

Model shifted_model(Model& params, std::vector<real_t> Delta);

double logcauchyv(double mode, double width);

double uniform(double mean, double width);

Model get_neighbour(Model& model, double w);

Model annealing_min(std::function<real_t(Model &model)> objective,
                    Model initial_guess);

Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess, int resolution);

Model brute_force_min(std::function<real_t(Model &model)> objective,
                      Model initial_guess);

Eigen::MatrixXd gradient_log(std::function<real_t(Model&)> objective,
                             Model point);

Eigen::MatrixXd compute_hessian(std::function<real_t(Model&)> objective,
                                Model point);

Model gradient_min(std::function<real_t(Model& model)> objective,
                   Model initial_guess);

Model mixed_min(std::function<real_t(Model& model)> objective,
                Model initial_guess);

Model mixed_min_8(std::function<real_t(Model& model)> objective,
                Model initial_guess);

Model mixed_min_16(std::function<real_t(Model& model)> objective,
                Model initial_guess);

Model skip_minimisation(std::function<real_t(Model& model)> objective,
                        Model initial_guess);

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

real_t save_data_compute_maximum(Epidata_t &all_times);

void draw_level_sets(std::function<real_t(Model &model)> objective,
                     Model point, int x_axis, int y_axis,
                     real_t x_range=10.0f, real_t y_range=10.0f,
                     size_t lines=16);

void draw_3dsurface(std::function<real_t(Model &model)> objective,
                    Model origin, int x_axis, int y_axis,
                    real_t x_range=10.0f, real_t y_range=10.0f,
                    size_t lines=16);

void sample_voxel_cube(float* buffer,
                       std::function<real_t(Model &model)> objective,
                       struct BoundingBox bounding_box,
                       real_t offset, real_t s2max,
                       size_t first_sample, size_t last_sample);

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

void guess_parameters_germline(Model &ground_truth, GuesserConfig options,
                               Model (*method_min)(std::function<real_t(Model&)>, Model));

Estimate get_estimate(real_t binwidth, size_t reference_pop,
                      Histogram_t incidence,
                      Model (*method_min)(std::function<real_t(Model&)>, Model));

void generate_histogram(Model &ground_truth, size_t seed, size_t dataset_size,
                        const std::vector<size_t> &end_nodes,
                        real_t &max_age, size_t &reference_pop,
                        real_t &binwidth, Histogram_t &incidence);

void guess_parameters(Model &ground_truth, GuesserConfig options,
                      Model (*method_min)(std::function<real_t(Model&)>, Model));

#endif
