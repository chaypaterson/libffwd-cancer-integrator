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
        std::vector<int> m_initial_pops;

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
            std::vector<int> metapopulation(m_stages, 0);
            m_initial_pops = metapopulation;
        }
};

Model product(const Model& graph_a, const Model& graph_b) {
    // compute and return the Kronecker product of the two models
    size_t n_vertices = graph_a.m_stages * graph_b.m_stages;
    Model graph_c(n_vertices);

    for (size_t vertex_a = 0; vertex_a < graph_a.m_stages; 
         ++vertex_a) {
        for (size_t vertex_b = 0; vertex_b < graph_b.m_stages; 
             ++vertex_b) {
            size_t vertex_c = vertex_a * graph_b.m_stages + vertex_b;
            // by default, assume no epistasis:
            graph_c.m_birth[vertex_c] = graph_a.m_birth[vertex_a]
                                            + graph_b.m_birth[vertex_b];
            graph_c.m_death[vertex_c] = graph_a.m_death[vertex_a]
                                            + graph_b.m_death[vertex_b];
            // FIXME I think this might be incorrect:
            graph_c.m_immig_rates[vertex_c] = graph_a.m_immig_rates[vertex_a]
                                            + graph_b.m_immig_rates[vertex_b];
            // TODO mutation rates here:
            // ... 
        }
    }

    return graph_c;
}

#include <iostream>
int main() {
    Model graph(3);

    std::cout << "(";
    for (auto& node : graph.m_birth)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;

    std::cout << "(";
    for (auto& node : graph.m_death)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;

    std::cout << "(";
    for (auto& node : graph.m_immig_rates)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;
        
    return 0;
}
