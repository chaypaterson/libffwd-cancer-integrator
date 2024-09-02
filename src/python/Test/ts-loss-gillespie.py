import pyffwd

def main(seed, total_runs):
    
    num_thr = max(1, os.cpu_count() - 2)
    runs_per_thr = runs // num_thr

    # System coefficients:
    rloh = 5e-7
    mu = 5e-8

    # Initialize the model
    model = pyffwd.Model(5)
    
    # Convert lists to RealVector
    birth_rates = pyffwd.list_to_vector([0, 0.05, 0.03, 0, 0])
    death_rates = pyffwd.list_to_vector([0, 0, 0, 0, 0])
    initial_pops = pyffwd.list_to_vector([1e6, 0, 0, 0, 0])
    
    # Set model parameters
    model.m_birth = birth_rates
    model.m_death = death_rates
    model.m_initial_pops = initial_pops

    # Define migration rates
    model.m_migr = [
        {1: mu, 2: rloh},         # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 4: 0.5 * rloh},  # From vertex 1 to 3 and 4
        {4: 0.5 * mu},            # From vertex 2 to 4
        {},                       # No migration from vertex 3
        {}                        # No migration from vertex 4
    ]

    # Define final vertices
    final_vertices = pyffwd.convert_to_vector_int([3, 4])
    
    # Initialize RNG and run simulations
    rng = pyffwd.GSL_RNG(seed)
    
    all_times = []

    def run_simulation(thread_id):
        rng_seed = seed + thread_id
        results = []
        times_to_final_vertices(model, rng_seed, runs_per_thr, final_vertices, results)
        return results

    # Run simulations in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_thr) as executor:
        futures = [executor.submit(run_simulation, i) for i in range(num_thr)]
        for future in concurrent.futures.as_completed(futures):
            all_times.extend(future.result())

    # Sorting results
    all_times.sort()

    # Print results for both types:
    for vertex in final_vertices:
        print(f"Type {vertex}:")

        mutant_times = [t for t, v in all_times if v == vertex]

        # Guard against invalid access:
        if len(mutant_times) < 1:
            print("No results")
            continue

        # Kaplan-Meier plot:
        time_max = max(t for t, _ in all_times)
        print("age, p1, p2,")
        print_kaplan_meier(time_max, mutant_times, len(all_times))

if __name__ == "__main__":
    import sys
    import os
    import concurrent.futures

    if len(sys.argv) < 3:
        print("Call this program with\n ./tsgillespie seed runs\n")
        sys.exit(1)

    seed = int(sys.argv[1])
    runs = int(sys.argv[2])

    main(seed, runs)
