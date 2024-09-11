import pyffwd
import sys

'''
Test to check the random numbers
def test_rng(seed, count=10):
    # Initialize RNG with seed
    rng = pyffwd.GSL_RNG(seed)
    # Generate random numbers using the appropriate method
    # Using 'uniform()' for [0, 1):
    numbers = [rng.uniform() for _ in range(count)]
    return numbers
'''

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

    #  set birth, death rates and initial_pops
    model = pyffwd.Model.create_model(
        5, 
        birth_rates=[0, 0.05, 0.03, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e6, 0, 0, 0, 0]
    )
    model.m_migr = [
        {1: mu, 2: rloh},         # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 4: 0.5 * rloh},  # From vertex 1 to 3 and 4
        {4: 0.5 * mu},            # From vertex 2 to 4
        {},                       # No migration from vertex 3
        {}                        # No migration from vertex 4
    ]
        
    final_vertices = pyffwd.list_vector_int([3, 4])
    
    '''
    # Test RNG consistency
    print("Python output with seed", seed, ":")
    rng_numbers = test_rng(seed)
    for number in rng_numbers:
        print(number)
    '''

    # Run some simulations and store the time and final node in all_times
    r = pyffwd.GSL_RNG(seed)
    all_times = []

    sim_state = pyffwd.GillespieInstance(model)

    # Make sure all_times is a list of tuples
    for _ in range(runs):
        time_result = pyffwd.first_passage_time_multiple(r, model,
                                                         final_vertices)
        all_times.append(time_result)

    # pass the list of tuples to times_to_final_vertices
    pyffwd.times_to_final_vertices(model, seed, runs, final_vertices, all_times)

    print("age, node,")
    for pair in all_times:
        age_formatted = f"{pair[0]:.3f}".rstrip('0').rstrip('.')
        node = pair[1]
        print(f"{age_formatted}, {node},")

    return 0

if __name__ == "__main__":
    main()
