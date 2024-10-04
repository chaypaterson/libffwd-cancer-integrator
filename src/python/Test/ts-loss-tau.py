import pyffwd
import numpy as np
import os
import threading

def times_to_final_vertices_tau(model, seed, runs_per_thr,
                            final_vertices, tau, result_container):
    rng = pyffwd.GSL_RNG(seed)
    for _ in range(runs_per_thr):
        result = pyffwd.first_passage_time_tau(
            rng, model, final_vertices, tau)
        result_container.append(result)

def print_kaplan_meier(time_max, mutant_times, sample_size):
    pyffwd.print_kaplan_meier(time_max,
                              pyffwd.RealVector(mutant_times),
                              sample_size)

if __name__ == "__main__":
      # default values
    sample_size = 10
    seed = 1
    tau = 0.1

    import sys
    if len(sys.argv) < 4:
        print("Call this program with\n /ts_tau seed runs tau \n")
        sys.exit(1)

    # TODO multithreading can break this test
    seed = int(sys.argv[1])
    sample_size = int(sys.argv[2])
    tau = float(sys.argv[3])

    # System coefficients:
    rloh = 5e-7
    mu = 5e-8
    print("%.4f" % rloh, "%.4f" % mu, "%.4f" % tau)

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
    times_to_final_vertices_tau(model, seed, sample_size,
                            final_vertices, tau, all_times)

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

        # Print the Kaplan-Meier
        print("age, p1, p2,", flush=True)
        print_kaplan_meier(time_max, mutant_times, sample_size)
        print()

        for time in mutant_times:
            print("%.8f" % time)

        print()
