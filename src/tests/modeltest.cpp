// Unit test for graph-model-spec.hpp:
#include <graph-model-spec.hpp>
#include <iostream>

int main() {
    using namespace clonal_expansion;

    Model graph(3);

    graph.m_birth[0] = 0.1;
    graph.m_birth[1] = 0.2;
    graph.m_birth[2] = 0.3;

    graph.m_death[0] = 0.1;
    graph.m_death[1] = 0.1;
    graph.m_death[2] = 0.05;

    graph.m_migr[0][1] = 1e-8;
    graph.m_migr[1][2] = 1e-7;

    graph.m_immig_rates[0] = 0.5;

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

    std::cout << "(";
    for (auto& node : graph.m_initial_pops)
        std::cout << node <<", ";
    std::cout << ")" << std::endl;

    std::cout << "Testing write out:" << std::endl;
    save_model(graph, "test.model");

    std::cout << "Testing read in:" << std::endl;
    Model graph2 = load_model("test.model");
    save_model(graph2, "test2.model");

    return 0;
}

