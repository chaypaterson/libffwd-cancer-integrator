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
    model.m_migr[0][1] = mu
    model.m_migr[0][2] = rloh
    model.m_migr[1][3] = 0.5 * mu
    model.m_migr[1][4] = 0.5 * rloh
    model.m_migr[2][4] = 0.5 * mu

    # Convert lists to RealVector
    model.m_birth = pyffwd.list_to_vector([0, 0.05, 0.03, 0, 0])
    model.m_death = pyffwd.list_to_vector([0, 0, 0, 0, 0])
    model.m_initial_pops = pyffwd.list_to_vector([1e6, 0, 0, 0, 0])

    final_vertices = [3, 4]

    # Run some simulations and store the time and final node in all_times:
    rng_wrapper = pyffwd.RNGWrapper(seed)
    all_times = rng_wrapper.first_passage_time_multiple(model, final_vertices)

    print("age, node,")
    for pair in all_times:
        print(f"{pair[0]}, {pair[1]},")

    return 0

if __name__ == "__main__":
    main()
