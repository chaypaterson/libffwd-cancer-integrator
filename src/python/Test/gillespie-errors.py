import pyffwd
import numpy as np
import math
import sys

def generate_dataset(seed, runs):
    # System coefficients
    rloh = 5e-7
    mu = 5e-8

    # Initialize the model with 5 compartments
    model = pyffwd.Model.create_model(
        5,
        birth_rates=[0, 0.05, 0.03, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e6, 0, 0, 0, 0]
    )

    # Set migration rates between compartments
    model.m_migr = [
        {1: mu, 2: rloh},  # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 
         4: 0.5 * rloh},  # From vertex 1 to 3 and 4
        {4: 0.5 * mu},    # From vertex 2 to 4
        {},                # No migration from vertex 3
        {}                 # No migration from vertex 4
    ]

    final_vertices = pyffwd.IntVector([3, 4])

    # Container for simulation results
    all_times = []

    # Set up random number generator
    rng = pyffwd.GSL_RNG(seed)

    # Run simulations using first_passage_time_multiple
    for _ in range(runs):
        result = pyffwd.first_passage_time_multiple(
            rng, model, final_vertices
        )
        all_times.append(result)

    # Flip results to map by final tumour type
    all_times_flipped = {3: [], 4: []}
    for time, vertex in all_times:
        all_times_flipped[vertex].append(time)

    return all_times_flipped


def compute_error(all_times_1, all_times_2, runs):
    # Compute maximum age
    age_max = 0
    for dataset in [all_times_1, all_times_2]:
        for times in dataset.values():
            if times:
                age_max = max(age_max, max(times))
            times.sort()

    # Reference populations
    reference_pop_1 = len(all_times_1[3]) + len(all_times_1[4])
    reference_pop_2 = len(all_times_2[3]) + len(all_times_2[4])

    # Parameters for sampling points
    num_sample_points = 256
    dt = age_max / num_sample_points

    # Output errors over time
    for age in np.arange(0, age_max + dt, dt):
        total_error = 0
        for tumour_type in [3, 4]:
            times_1 = (
                pyffwd.RealVector(all_times_1[tumour_type]) 
                if all_times_1[tumour_type] 
                else pyffwd.RealVector()
            )
            times_2 = (
                pyffwd.RealVector(all_times_2[tumour_type]) 
                if all_times_2[tumour_type] 
                else pyffwd.RealVector()
            )

            s1 = (1.0 if not times_1 
                  else pyffwd.surv_kaplan_meier(
                      age, times_1, reference_pop_1))
            s2 = (1.0 if not times_2 
                  else pyffwd.surv_kaplan_meier(
                      age, times_2, reference_pop_2))
            
            error = (s1 - s2) ** 2
            total_error += error

        total_error = math.sqrt(total_error / 2)
        print(f"{age:.20f},{total_error:.20f},")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Call this program with: ")
        print("    python gillespie_sampler.py seed runs")
        sys.exit(1)

    seed = int(sys.argv[1])
    runs = int(sys.argv[2])

    # Generate two datasets with different seeds
    all_times_1 = generate_dataset(seed, runs)
    all_times_2 = generate_dataset(seed + runs, runs)

    # Compute and output the error between the two datasets
    compute_error(all_times_1, all_times_2, runs)
