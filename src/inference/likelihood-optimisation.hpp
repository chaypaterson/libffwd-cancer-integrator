#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <cstdlib>
#include <cstring>

enum Minimise {ANNEALING, GRADIENT};

namespace clonal_expansion {

class GuesserConfig 
{
public:
    // FLAG AND ARGUMENT OPTIONS:
    size_t seed;
    size_t dataset_size;
    // TODO load histogram from file
    // SINGLE FLAG OPTIONS:
    bool include_germline = false; // mixed germline/sporadic study or not (default not)
    bool resample_after   = false; // resample or not (default not)
    // what minimisation method to use (default annealing, can do gradient)
    enum Minimise minimise_with = ANNEALING; // TODO could just be a bool
    bool level_sets       = false; // whether or not to draw with level sets
    bool draw_mesh        = false; // whether or not to draw 3d plots of the likelihood

    GuesserConfig(int argc, char* argv[])
    {
        // for arguments that come in pairs: " --seed 5 " etc.
        for (int arg = 1; arg < argc - 1; ++arg) {
            if (!strcmp(argv[arg], "--seed"))
                seed = atoi(argv[arg + 1]);
            if (!strcmp(argv[arg], "--sample_size"))
                dataset_size = atoi(argv[arg + 1]);
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
