#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <cstdlib>
#include <cstring>

enum Minimise {ANNEALING, GRADIENT};

namespace clonal_expansion {

/* GuesserConfig interprets command-line arguments as configuration options.
 * To add a command-line argument:
 * 1. add a public member variable to GuesserConfig
 * 2. make sure the member variable has a default value
 * 3. add the matching "if (match) variable = ...;" to the constructor below
 *
 * TODO a better way of iterating over the variable-string pairs would be
 * good
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
    enum Minimise minimise_with = ANNEALING; // TODO could just be a bool
    bool level_sets       = false; // whether or not to draw with level sets
    bool draw_mesh        = false; // whether or not to draw 3d plots of the likelihood

    // The constructor:
    inline GuesserConfig(int argc, char* argv[])
    {
        // for arguments that come in pairs: " --seed 5 " etc.
        for (int arg = 1; arg < argc - 1; ++arg) {
            if (!strcmp(argv[arg], "--seed"))
                seed = atoi(argv[arg + 1]);
            if (!strcmp(argv[arg], "--sample_size"))
                dataset_size = atoi(argv[arg + 1]);
            if (!strcmp(argv[arg], "--mesh_lines"))
                mesh_lines = atoi(argv[arg + 1]);
            if (!strcmp(argv[arg], "--mesh_x_range")) 
                mesh_x_range = atof(argv[arg + 1]);
            if (!strcmp(argv[arg], "--mesh_y_range")) 
                mesh_y_range = atof(argv[arg + 1]);
        }
        // for arguments that are just isolated flags: "--annealing" etc.
        for (int arg = 1; arg < argc; ++arg) {
            if (!strcmp(argv[arg], "--with_germline")) include_germline = true;
            if (!strcmp(argv[arg], "--resample")) resample_after = true;
            if (!strcmp(argv[arg], "--annealing")) minimise_with = ANNEALING;
            if (!strcmp(argv[arg], "--gradient")) minimise_with = GRADIENT;
            if (!strcmp(argv[arg], "--draw_level_sets")) level_sets = true;
            if (!strcmp(argv[arg], "--draw_3d_meshes")) draw_mesh = true;
        }
    }
private:
};

} // namespace clonal_expansion

#endif //LIKELIHOOD_H
