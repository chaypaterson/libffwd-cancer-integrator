#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <cstdlib>
#include <cstring>

enum Minimise {ANNEALING, GRADIENT};

class GuesserConfig 
{
public:
    size_t seed;
    size_t dataset_size;
    bool include_germline = false; // mixed germline/sporadic study or not (default not)
    bool resample_after   = false; // resample or not (default not)
    enum Minimise minimise_with = ANNEALING;
    // what minimisation method to use (default annealing, can do gradient)

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
        }
    }
private:
};

#endif //LIKELIHOOD_H
