#include <vector>
#include <map>
#include <cstddef>
#ifndef GRAPH_MODEL_DEF
#define GRAPH_MODEL_DEF

/* A class for containing a specification of a birth-death-mutation-immigration
 * process on a graph.
 * This class contains the parameters and initial conditions of a model.
 * The parameters are the birth rates, death rates, and immigration rates for
 * each vertex, and the migration/mutation rate graph. The migration graph is
 * represented as a vector of maps for efficiency, but can equally be thought of
 * as a matrix.
 *
 * This class should be agnostic to the method used to solve the model: it
 * should be possible to feed instances of Model into Gillespie algorithm
 * solvers, as well as new methods.
 */

typedef double real_t;

class Model {
    public:
        size_t m_stages{0};
        std::vector<real_t> m_birth;
        std::vector<real_t> m_death;
        std::vector<real_t> m_immig_rates;
        std::vector<std::map<int, real_t>> m_migr;
        std::vector<real_t> m_initial_pops;

        // constructor:
        Model(size_t n_vertices) {
            m_stages = n_vertices;
            std::vector<real_t> nodes(m_stages, 0);
            m_birth = nodes;
            m_death = nodes;
            m_immig_rates = nodes;
            // Initialise the graph of mutation rates:
            std::map<int, real_t> neighbours;
            std::vector<std::map<int, real_t>> edges(m_stages, neighbours);
            m_migr = edges;
            std::vector<real_t> metapopulation(m_stages, 0);
            m_initial_pops = metapopulation;
        }
};

// end of header guard GRAPH_MODEL_DEF
#endif
