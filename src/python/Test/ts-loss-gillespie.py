import pyffwd
import numpy as np
import os
import threading

def times_to_final_vertices(model, seed, runs_per_thr,
                            final_vertices, result_container):
    rng = pyffwd.GSL_RNG(seed)
    for _ in range(runs_per_thr):
        result = pyffwd.first_passage_time_multiple(
            rng, model, final_vertices)
        result_container.append(result)

def print_kaplan_meier(time_max, mutant_times, study_population):
    pyffwd.print_kaplan_meier(time_max,
                              pyffwd.RealVector(mutant_times),
                              study_population)

if __name__ == "__main__":
      # default values
    runs_per_thr = 10
    seed = 1

    import sys
    if len(sys.argv) < 3:
        print("Call this program with\n ./tsgillespie seed runs \n")
        sys.exit(1)

    seed = int(sys.argv[1])
    total_runs = int(sys.argv[2])
    
    num_thr = max(1, os.cpu_count() - 2)  
    runs_per_thr = total_runs // num_thr

    # System coefficients:
    rloh = 5e-7
    mu = 5e-8

    # Define the model
    model = pyffwd.Model.create_model(
        5, 
        birth_rates=[0, 0.05, 0.03, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e6, 0, 0, 0, 0]
    )
    model.m_migr = [
        {1: mu, 2: rloh},
        {3: 0.5 * mu, 4: 0.5 * rloh},
        {4: 0.5 * mu},
        {},
        {}
    ]

    final_vertices = pyffwd.IntVector([3, 4])

    all_times = []

    times = [[] for _ in range(num_thr)]
    threads = []

    for i in range(num_thr):
        thread = threading.Thread(
            target=times_to_final_vertices, 
            args=(model, seed + i, runs_per_thr,
                  final_vertices, times[i])
        )
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    for thread_results in times:
        all_times.extend(thread_results)

    all_times.sort(key=lambda x: x[0])

    study_population = len(all_times)
    time_max = all_times[-1][0] if all_times else 0

    # Print Kaplan-Meier
    for tumour_type in final_vertices:
        print(f"Type {tumour_type}:")

        mutant_times = [result[0] for result in all_times 
                        if result[1] == tumour_type]

        if not mutant_times:
            print("No results")
            sys.exit(2)

        # Print the Kaplan-Meier curve
        print("age, p1, p2,", flush=True)
        print_kaplan_meier(time_max, mutant_times, study_population)
        print()
