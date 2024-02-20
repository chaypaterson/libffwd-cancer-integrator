#ifndef MOST_LIKELY_DEF
#define MOST_LIKELY_DEF
#include "graph-model-spec.hpp"
#include "fast-forward.hpp"
#include <vector>
#include <cmath>
#include <cstdio>

/* Functions helpful for computing the likelihood function for branching
 * processes on graphs.
 */

// TODO replace real_t with double

// A list of ages and (integer) cancer subtypes:
typedef std::vector<std::pair<gmsce::real_t, int>> epidata_t;

// A function pointer type with a specific signature:
typedef gmsce::real_t mappable_t(gmsce::real_t age, int node, gmsce::real_t prob, gmsce::real_t dprob);

// Map a function onto the ages and nodes of the dataset, as well as the
// computed probability and hazard:
void map_onto_data(gmsce::Model& params, const epidata_t& this_data, 
                   mappable_t *mapme, gmsce::real_t *result);

// Functions we want to map for our tests:
mappable_t print_test;
mappable_t logdprob;

#endif // MOST_LIKELY_DEF
