import numpy as np
import multiprocessing as mp
import pyffwd

def simulate_times(seed, runs_per_thread, final_vertices_list, rloh, mu):
    # Initialize the model inside the worker function
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

    # Convert the list to a pyffwd.IntVector
    final_vertices = pyffwd.IntVector(final_vertices_list)
    
    # Initialize an empty pyffwd.RealVector
    times = pyffwd.RealVector()
    
    # Call the simulation function
    pyffwd.times_to_final_vertices_poly(model, seed, runs_per_thread, final_vertices, times)
    
    # Convert the RealVector to a Python list
    return list(times)

def main():
    # Default values
    seed = 1
    runs_per_thr = int(1e7)  # Number of runs per thread
    
    # If not enough command line arguments, display the usage message
    import sys
    if len(sys.argv) < 3:
        print("Usage: python ts_gillespie.py <seed> <runs>")
        sys.exit(1)

    # Parse arguments from the command line
    seed = int(sys.argv[1])
    total_runs = int(sys.argv[2])
    
    # Number of threads (set to available CPUs minus 2)
    num_threads = max(1, mp.cpu_count() - 2)
    runs_per_thr = total_runs // num_threads

    # System coefficients
    rloh = 5e-7
    mu = 5e-8

    final_vertices_list = [3, 4]

    # Run simulations across multiple threads
    all_times = []
    with mp.Pool(num_threads) as pool:
        # Dispatch simulation jobs to the thread pool
        results = [pool.apply_async(simulate_times, 
                                    args=(seed + i, runs_per_thr, final_vertices_list, rloh, mu)) 
                   for i in range(num_threads)]
        
        # Collect results from all threads
        for result in results:
            all_times.extend(result.get())

    # Sort all simulation times and determine the total study population
    all_times.sort()
    study_population = len(all_times)
    time_max = all_times[-1]

    all_times_real_vector = pyffwd.RealVector(all_times)
    
    print("age, S12, p12,", flush=True)
    # Print the Kaplan-Meier plot
    pyffwd.print_kaplan_meier(time_max, all_times_real_vector, study_population)
    print()

if __name__ == '__main__':
    main()
