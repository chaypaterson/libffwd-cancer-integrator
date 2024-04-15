#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <cstdlib>
#include <cstring>

namespace clonal_expansion {

/* GuesserConfig interprets command-line arguments as configuration options.
 * To add a command-line argument:
 * 1. add a public member variable to GuesserConfig
 * 2. make sure the member variable has a default value
 * 3. add the matching setter:
 *       set_...(cmdarg, "--...", ...); 
 *    to the constructor below
 */

class GuesserConfig 
{
public:
    // FLAG AND ARGUMENT OPTIONS (paired options):
    size_t seed = 1; // RNG seed
    size_t dataset_size = 10; // sample size
    size_t mesh_lines = 16; // resolution of the heatmap plot of the likelihood
    // the heatmap plot will run from best_guess/these to best_guess * these:
    double mesh_x_range = 10.0f; 
    double mesh_y_range = 10.0f;
    // TODO load histogram from file
    // SINGLE FLAG OPTIONS (boolean options):
    bool include_germline = false; // mixed germline/sporadic study or not (default not)
    bool resample_after   = false; // resample or not (default not)
    // what minimisation method to use (default annealing, can do gradient)
    bool minimise_with_gradient = false;
    bool level_sets       = false; // whether or not to draw with level sets
    bool draw_mesh        = false; // whether or not to draw 3d plots of the likelihood

    // The constructor:
    inline GuesserConfig(int argc, char* argv[])
    {
        // for arguments that come in pairs: " --seed 5 " etc.
        for (char* *cmdarg = argv; cmdarg < argv + argc - 1; ++cmdarg) {
            set_pair(cmdarg, "--seed", seed);
            set_pair(cmdarg, "--sample_size", dataset_size);
            set_pair(cmdarg, "--mesh_lines", mesh_lines);
            set_pair(cmdarg, "--mesh_x_range", mesh_x_range);
            set_pair(cmdarg, "--mesh_y_range", mesh_y_range);
        }
        // for arguments that are just isolated flags: "--annealing" etc.
        for (char* *cmdarg = argv; cmdarg < argv + argc; ++cmdarg) {
            set_bool(cmdarg, "--with_germline", include_germline);
            set_bool(cmdarg, "--with_germline", include_germline);
            set_bool(cmdarg, "--resample", resample_after);
            set_bool(cmdarg, "--gradient", minimise_with_gradient);
            set_bool(cmdarg, "--draw_level_sets", level_sets);
            set_bool(cmdarg, "--draw_3d_meshes", draw_mesh);
        }
    }
private:
    inline void set_pair(char* *cmdarg, const char key[], size_t &member);
    inline void set_pair(char* *cmdarg, const char key[], double &member);
    inline void set_bool(char* *cmdarg, const char key[], bool &member);
};

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
