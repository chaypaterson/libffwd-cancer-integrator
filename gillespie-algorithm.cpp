#include <graph-model-spec.hpp>

// TODO a class defining an iterator for the Gillespie algorithm 
// this should contain:
// data:
//    * a time variable
//    * a vector of populations
//    * a Model specification (the rate parameters)
// methods:
//    * a method to compute Gamma
//    * a method to accept a variable x < Gamma and return an event, which
//      consists of a change in the populations
//    * methods to generate random variables x and Deltat
//    * a method to perform one step, changing the populations and advancing time

class Gillespie {
    public:
        double m_time{0};
        size_t m_n_populations{0};
        std::vector<int> m_metapopulation;
        Model m_parameters;

        Gillespie(Model parameters) {
            m_parameters = parameters;
            m_n_vertices = parameters.m_n_vertices;
            m_metapopulation = parameters.m_initial_popns;
        }
        
        //void step() {
        //}
    private:
        // rng;
        // double Gamma();
        // double x();
        // double Deltat();
        // std::vector<int> delta_popns();
};
