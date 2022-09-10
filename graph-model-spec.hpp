#include <vector>
#include <map>
// A class for containing a specification of a birth-death-mutation-immigration
// process on a graph.
// This class contains the parameters and initial conditions of a model.

class Model {
    public:
        size_t m_stages{0};
        std::vector<double> m_birth;
        std::vector<double> m_death;
        std::vector<double> m_immig_rates;
        std::vector<std::map<int, double>> m_migr;
        std::vector<double> m_initial_pops;

        // constructor:
        Model(size_t n_vertices) {
            m_stages = n_vertices;
            std::vector<double> nodes(m_stages, 0);
            m_birth = nodes;
            m_death = nodes;
            m_immig_rates = nodes;
            std::map<int, double> neighbours;
            std::vector<std::map<int, double>> edges(m_stages, neighbours);
            m_migr = edges;
            std::vector<double> metapopulation(m_stages, 0);
            m_initial_pops = metapopulation;
        }
};
