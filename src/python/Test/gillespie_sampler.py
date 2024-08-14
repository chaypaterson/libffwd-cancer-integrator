import pyffwd
import sys

def main():
    if len(sys.argv) < 3:
        print("Call this program with:")
        print("    python gillespie_sampler.py seed runs")
        return 1
    
    seed = int(sys.argv[1])
    runs = int(sys.argv[2])

    # System coefficients:
    rloh = 5.26337e-7
    mu = 2.16427e-8

    model = pyffwd.Model(5)
    model.m_migr = [
        {1: mu, 2: rloh}, # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 4: 0.5 * rloh}, # From vertex 1 to 3 and 4
        {4: 0.5 * mu} # From vertex 2 to 4
    ]

    # Convert lists to RealVector
    model.m_birth = pyffwd.list_to_vector([0, 0.05, 0.03, 0, 0])
    model.m_death = pyffwd.list_to_vector([0, 0, 0, 0, 0])
    model.m_initial_pops = pyffwd.list_to_vector([1e6, 0, 0, 0, 0])

    final_vertices = [3, 4]

    # Run some simulations and store the time and final node in all_times:
    r = pyffwd.seed_gsl_rng(seed)
    all_times = []

    sim_state = pyffwd.GillespieInstance(model) # causes segfault FIXME
    print("hello")
    all_times = pyffwd.first_passage_time_multiple(r, model, final_vertices)

    print("wololo")

    pyffwd.times_to_final_vertices(model, seed, runs, final_vertices, all_times)

    print("age, node,")
    for pair in all_times:
        print(f"{pair[0]}, {pair[1]},")

    return 0

if __name__ == "__main__":
    main()
