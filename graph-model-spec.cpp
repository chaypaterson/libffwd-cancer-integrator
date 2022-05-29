#include <vector>
#include <map>
// A class for containing a specification of a birth-death-mutation-immigration
// process on a graph.
// This class contains the parameters and initial conditions of a model.

class Model {
    public:
        size_t m_n_vertices{0};
        std::vector<double> m_birth_rates;
        std::vector<double> m_death_rates;
        std::vector<double> m_immig_rates;
        std::vector<std::map<int, double>> m_mutations;

        // constructor:
        Model(size_t n_vertices) {
            m_n_vertices = n_vertices;
            std::vector<double> nodes(m_n_vertices, 0);
            m_birth_rates = nodes;
            m_death_rates = nodes;
            m_immig_rates = nodes;
            std::map<int, double> neighbours;
            std::vector<std::map<int, double>> edges(m_n_vertices, neighbours);
            m_mutations = edges;
        }
};

Model product(const Model& graph_a, const Model& graph_b) {
    // compute and return the Kronecker product of the two models
    size_t n_vertices = graph_a.m_n_vertices * graph_b.m_n_vertices;
    Model graph_c(n_vertices);

    for (size_t vertex_a = 0; vertex_a < graph_a.m_n_vertices; 
         ++vertex_a) {
        for (size_t vertex_b = 0; vertex_b < graph_b.m_n_vertices; 
             ++vertex_b) {
            size_t vertex_c = vertex_a * graph_b.m_n_vertices + vertex_b;
            // by default, assume no epistasis:
            graph_c.m_birth_rates[vertex_c] = graph_a.m_birth_rates[vertex_a]
                                            + graph_b.m_birth_rates[vertex_b];
            graph_c.m_death_rates[vertex_c] = graph_a.m_death_rates[vertex_a]
                                            + graph_b.m_death_rates[vertex_b];
            // FIXME I think this might be incorrect:
            graph_c.m_immig_rates[vertex_c] = graph_a.m_immig_rates[vertex_a]
                                            + graph_b.m_immig_rates[vertex_b];
        }
    }

    return graph_c;
}

#include <iostream>
int main() {
    Model graph(3);

    std::cout << "(";
    for (auto& node : graph.m_birth_rates)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;

    std::cout << "(";
    for (auto& node : graph.m_death_rates)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;

    std::cout << "(";
    for (auto& node : graph.m_immig_rates)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;
        
    return 0;
}
