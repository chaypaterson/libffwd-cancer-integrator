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
