import pyffwd

def main(seed, total_runs):
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
    model.m_migr[0][1] = 5e-8
    model.m_migr[0][2] = 5e-7
    model.m_migr[1][3] = 0.5 * 5e-8
    model.m_migr[1][4] = 0.5 * 5e-7
    model.m_migr[2][4] = 0.5 * 5e-8

    # Define final vertices
    final_vertices = [3, 4]
    
    # Initialize RNG and run simulations
    rng = pyffwd.RNGWrapper(seed)
    
    # Debug: Print parameters
    print(f"Running simulations with seed={seed}, final_vertices={final_vertices}, total_runs={total_runs}")
    
    all_times = []
    for _ in range(total_runs):
        times = rng.first_passage_time_multiple(model, final_vertices)
        if times:
            all_times.extend(times)
        else:
            print("No results returned from first_passage_time_multiple for one of the runs.")
    
    # Debug: Check the results
    if not all_times:
        print("No results returned from all runs.")
        return
    
    print(f"Number of results: {len(all_times)}")

    # Sort and process results
    all_times.sort()
    if not all_times:
        print("No results to process after sorting.")
        return
    
    time_max = max([t[0] for t in all_times])

    # Print results for each type
    for type in final_vertices:
        print(f"Type {type}:")
        mutant_times = [t[0] for t in all_times if t[1] == type]
        
        if not mutant_times:
            print("No results for this type")
            continue
        
        # Convert mutant_times to RealVector
        mutant_times_vector = pyffwd.list_to_vector(mutant_times)
        
        # Print Kaplan-Meier estimate
        print("age, p1, p2,")
        pyffwd.print_kaplan_meier(time_max, mutant_times_vector, len(all_times))

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python ts-loss-gillespie.py <seed> <runs>")
        sys.exit(1)
    
    seed = int(sys.argv[1])
    total_runs = int(sys.argv[2])
    main(seed, total_runs)
