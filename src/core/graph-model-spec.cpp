#include <vector>
#include <map>

#include <graph-model-spec.hpp>
// A class for containing a specification of a birth-death-mutation-immigration
// process on a graph.
// This file contains a unit test.

namespace clonal_expansion {

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

// serialisation and deserialisation:
void save_model(Model graph, std::string filename) {
    std::ofstream savefile(filename, std::ios::out);
    // File header: magic numbers
    savefile << "GMSCE" << std::endl;
    // Number of nodes, needed for the constructor:
    savefile << graph.m_stages << std::endl;

    // On the next line, write the birth rates:
    for (auto& elem : graph.m_birth) savefile << elem << " ";
    savefile << std::endl;
    // death rates:
    for (auto& elem : graph.m_death) savefile << elem << " ";
    savefile << std::endl;
    // and immigration rates:
    for (auto& elem : graph.m_immig_rates) savefile << elem << " ";
    savefile << std::endl;
    // Put the initial pops here for conceptual consistency:
    for (auto& elem : graph.m_initial_pops) savefile << elem << " ";
    savefile << std::endl;

    // The migration matrix is the rest of the file, and varies in size:
    size_t in_vertex = 0;
    for (auto& line : graph.m_migr) {
        for (auto const& pair : line) {
            savefile << in_vertex << " ";
            savefile << pair.first << " ";
            savefile << pair.second << std::endl;
        }
        ++in_vertex;
    }
}

Model load_model(std::string filename) {
    std::ifstream savefile(filename, std::ios::in);
    std::string header;
    savefile >> header;
    // Check header is correct:
    // TODO test
    std::string testheader{"GMSCE"};
    if (testheader.compare(header)) {
        std::cout << "Invalid model file header: " << header;
        std::cout << std::endl;
        exit(1);
    }

    size_t n_vertices;
    savefile >> n_vertices;
    // Construct model with all parameters zero initialised:
    Model newmodel(n_vertices);

    // Set vectors:
    for (auto& elem : newmodel.m_birth) savefile >> elem;
    for (auto& elem : newmodel.m_death) savefile >> elem;
    for (auto& elem : newmodel.m_immig_rates) savefile >> elem;
    for (auto& elem : newmodel.m_initial_pops) savefile >> elem;

    // Set graph:
    while (!savefile.eof()) {
        size_t in_vertex, out_vertex;
        real_t rate;
        savefile >> in_vertex >> out_vertex >> rate;
        newmodel.m_migr[in_vertex][out_vertex] = rate;
    }

    return newmodel;
}

}
