#include <stdio.h>
#define __USE_MISC
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "annealing.h"

/* annealing.c
 * This program implements a simple simulated annealing algorithm. It
 * demonstrates the use of simulated annealing to find the global minimum of a
 * test function with a shallow global minimum at e=2.71828...
 *
 * Chay Paterson (chay.paterson@manchester.ac.uk)
 * & Miaomiao Gao (miaomiao.gao@postgrad.manchester.ac.uk)
 */

double energy(double x) {
    double power = log(x);
    power /= x;
    return -2 * exp(power);
}

double uniform(double x, double w) {
    /* return a uniformly distributed random variable in the interval [x-w, x+w]
     */
    double rand_num = (double) rand() / RAND_MAX; 
    double scaled_num = (rand_num * 2 - 1) * w; 
    return x + scaled_num;
}

double normalv(double mean, double var) {
    /* return a normally distributed random variable centred on "mean" with
     * variance "var" */
    double rand_num = (double) rand() / RAND_MAX;
    double quantile = 1.0; /* TODO inverse erf? */
    return mean + quantile;
}

double cauchyv(double mean, double width) {
    /* return a Cauchy distributed random variable centred on "mean" with scale
     * coefficient "width" */
    double rand_num = (double)rand() / RAND_MAX;
    double quantile = width * tan(M_PI * (rand_num - 0.5f));
    return mean + quantile;
}

double annealing_min(double (*obj_fun)(double),
                     double x_initial_guess,
                     double precision)  {
    // Constants:
    const double Tmin = 1e-8; // Minimal temperature
    const double delta = 0.98; // The rate of temperature drop
    const double min_width = 1.0f;

    srand((unsigned)time(NULL));

    // Initialise variables:
    double Temp = 100;             // Initial temperature 
    double x = x_initial_guess;
    double best_x = x ;
    double best_y = obj_fun(best_x);
    double delta_x = 0;
    double old_best_y = best_y; // Used for adaptive width size

    // Width of neighbourhood:
    double w = Temp; // initial width

    // Main loop:
    while ((Temp > Tmin) || (delta_x * delta_x > precision * precision)) {
        double x_new = cauchyv(x, w);
        double y_new = obj_fun(x_new);

        // Don't count -NaN values as minima:
        if (isnan(y_new)) continue;

        double delta_y = y_new - best_y;

        if ((delta_y < 0) || ((rand() / (RAND_MAX + 1.0)) < exp(-delta_y / Temp))) {
            delta_x = x_new - best_x;
            x = x_new;
            best_x = x;
            best_y = y_new;
            printf("%g\n", x_new);
        }
        Temp *= delta;
        /* If best_y has improved by more than the current temperature, reduce
         * the width:*/
        if ((old_best_y - best_y) > Temp) {
            w = 0.5 * (min_width + w); // This will smoothly approach min_width
            old_best_y = best_y;
            printf("Shrinking width\n");
        }
    }
    return best_x;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("This program should be invoked like this:\n");
        printf("    ./anneal [precision]\n");
        printf("where [precision] is your chosen precision for the answer.\n");
        return 1;
    }

    double x_guess = 100.0;
    double precision = atof(argv[1]);
    printf("Minimising objective function... \n ");
    double x_min = annealing_min(energy, x_guess, precision);
    printf("%.17g %.17g\n", x_min, precision);
    printf("Done\n");

    return 0;
}

