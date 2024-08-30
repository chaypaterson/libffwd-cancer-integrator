#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <cstdlib>
#include <cstring>
#include <string>

namespace clonal_expansion {

/* GuesserConfig interprets command-line arguments as configuration options.
 * To add a command-line argument:
 * 1. add a public member variable to GuesserConfig
 * 2. make sure the member variable has a default value (unless it's a string)
 * 3. add the matching setter:
 *       set_...(cmdarg, "--...", ...); 
 *    to the constructor below
 */

/* TODO use inheritance to make the basic/private methods agnostic to the
   program? */
class GuesserConfig 
{
private:
    // factory for string arguments:
    std::string to_string(char* *argv, const char key[]);
public:
    // FLAG AND ARGUMENT OPTIONS (paired options):
    // integer options:
    size_t seed = 1; // RNG seed
    size_t dataset_size = 10; // sample size
    size_t mesh_lines = 16; // resolution of the heatmap plot of the likelihood
    size_t num_child_threads = 0; // parallelisation in jackknife_and_save
    size_t voxel_res = 16; // side resolution of 3d cube
    // float options:
    // the heatmap plot will run from best_guess/these to best_guess * these:
    double mesh_x_range = 10.0f; 
    double mesh_y_range = 10.0f;

    // string options:
    // NB the defaults for all string options are empty (""), this is set in the
    // to_string factory below.
    std::string histogram_file; // load histogram from file
    std::string estimate_file; // save (serialise) best guess to this file
    std::string voxel_file; // save voxel sampled function to this file

    // SINGLE FLAG OPTIONS (boolean options):
    bool include_germline = false; // mixed germline/sporadic study or not (default not)
    bool resample_after   = false; // resample or not (default not)
    // what minimisation method to use (default annealing, can do gradient):
    bool minimise_with_gradient = false;
    bool minimise_brute_force = false;
    bool minimise_with_mixed = false;
    bool minimise_with_mixed_8 = false;
    bool minimise_with_mixed_16 = false;
    bool skip_minimisation = false; // skip minimisation entirely
    bool level_sets       = false; // whether or not to draw with level sets
    bool draw_mesh        = false; // whether or not to draw 3d plots of the likelihood
    bool draw_voxels      = false; // draw voxels

    // The constructor:
    inline GuesserConfig(int argc, char* argv[]) :
        histogram_file(to_string(argv, "--load_histogram")),
        estimate_file(to_string(argv, "--save_estimate")),
        voxel_file(to_string(argv, "--voxel_file"))
    {
        // for arguments that come in pairs: " --seed 5 " etc.
        for (char* *cmdarg = argv; cmdarg < argv + argc - 1; ++cmdarg) {
            set_pair(cmdarg, "--seed", seed);
            set_pair(cmdarg, "--sample_size", dataset_size);
            set_pair(cmdarg, "--child_threads", num_child_threads);
            set_pair(cmdarg, "--mesh_lines", mesh_lines);
            set_pair(cmdarg, "--mesh_x_range", mesh_x_range);
            set_pair(cmdarg, "--mesh_y_range", mesh_y_range);
            set_pair(cmdarg, "--voxel_resolution", voxel_res);
        }
        // for arguments that are just isolated flags: "--annealing" etc.
        for (char* *cmdarg = argv; cmdarg < argv + argc; ++cmdarg) {
            set_bool(cmdarg, "--with_germline", include_germline);
            set_bool(cmdarg, "--resample", resample_after);
            set_bool(cmdarg, "--gradient", minimise_with_gradient);
            set_bool(cmdarg, "--with_brute_force", minimise_brute_force);
            set_bool(cmdarg, "--mixed_min", minimise_with_mixed);
            set_bool(cmdarg, "--mixed_8", minimise_with_mixed_8);
            set_bool(cmdarg, "--mixed_16", minimise_with_mixed_16);
            set_bool(cmdarg, "--skip_minimisation", skip_minimisation);
            set_bool(cmdarg, "--draw_level_sets", level_sets);
            set_bool(cmdarg, "--draw_3d_meshes", draw_mesh);
        }
    }
private:
    inline void set_pair(char* *cmdarg, const char key[], size_t &member);
    inline void set_pair(char* *cmdarg, const char key[], double &member);
    inline void set_bool(char* *cmdarg, const char key[], bool &member);
};

// TODO argv is null-terminated, so we don't need to pass argc
inline std::string GuesserConfig::to_string(char* *argv, const char key[]) {
    for (char* *cmdarg = argv; *(cmdarg + 1); // while next is not null...
         ++cmdarg) {
        if (!strcmp(*cmdarg, key)) {
            std::string outstring(cmdarg[1]);
            return outstring;
        }
    }
    // else:
    return ""; // empty string
}

inline void GuesserConfig::set_pair(char* *cmdarg, const char key[], size_t &member) {
    if (!strcmp(*cmdarg, key)) member = atoi(cmdarg[1]);
}

inline void GuesserConfig::set_pair(char* *cmdarg, const char key[], double &member) {
    if (!strcmp(*cmdarg, key)) member = atof(cmdarg[1]);
}

inline void GuesserConfig::set_bool(char* *cmdarg, const char key[], bool &member) {
    if (!strcmp(*cmdarg, key)) member = true;
}

} // namespace clonal_expansion

#endif //LIKELIHOOD_H
