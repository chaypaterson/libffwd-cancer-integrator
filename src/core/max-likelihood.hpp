#ifndef MOST_LIKELY_DEF
#define MOST_LIKELY_DEF
#include "graph-model-spec.hpp"
#include "flying-conjugates.hpp"
#include <vector>
#include <cmath>
#include <cstdio>

/* Functions helpful for computing the likelihood function for branching
 * processes on graphs.
 */

// A list of ages and (integer) cancer subtypes:
typedef std::vector<std::pair<real_t, int>> epidata_t;

// Map a function onto the ages and nodes of the dataset, as well as the
// computed probability and hazard:
typedef real_t mappable_t(real_t age, int node, real_t prob, real_t dprob);
void map_onto_data(Model& params, const epidata_t& this_data, 
                   mappable_t *mapme, real_t *result);

// Functions we want to map for our tests:
mappable_t print_test;
mappable_t logdprob;

// Unit tests:
void unit_test(Model& params, const epidata_t& all_times);
real_t loglikelihood_test(Model& params, const epidata_t& all_times);

#endif // MOST_LIKELY_DEF
