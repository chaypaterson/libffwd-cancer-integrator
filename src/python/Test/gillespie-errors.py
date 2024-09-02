import pyffwd
import numpy as np
import sys

def generate_dataset(seed, runs):
    # System coefficients:
    rloh = 5.0e-7
    mu = 5.0e-8

    model = pyffwd.Model(5)
    model.m_migr = [
        {1: mu, 2: rloh},
        {3: 0.5 * mu, 4: 0.5 * rloh},
        {4: 0.5 * mu},
        {},
        {}
    ]

    # Convert lists to RealVector
    model.m_birth = pyffwd.list_to_vector([0, 0.05, 0.03, 0, 0])
    model.m_death = pyffwd.list_to_vector([0, 0, 0, 0, 0])
    model.m_initial_pops = pyffwd.list_to_vector([1e6, 0, 0, 0, 0])
    
    final_vertices = pyffwd.convert_to_vector_int([3, 4])

    # Run some simulations and store the time and final node in all_times:
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

    all_times_flipped = {}
    for time, vertex in all_times:
        if vertex not in all_times_flipped:
            all_times_flipped[vertex] = []
        all_times_flipped[vertex].append(time)

    return all_times_flipped

def main(seed, runs):
    all_times_1 = generate_dataset(seed, runs)
    all_times_2 = generate_dataset(seed + runs, runs)

    age_max = 0
    for pair in all_times_1.values():
        if pair:  # Check if the list is not empty
            age_max = max(age_max, max(pair))
            pair.sort()
    
    for pair in all_times_2.values():
        if pair:
            age_max = max(age_max, max(pair))
            pair.sort()

    # Compute scaled root mean square difference in survival curves
    reference_pop_1 = len(all_times_1.get(3, [])) + len(all_times_1.get(4, []))
    reference_pop_2 = len(all_times_2.get(3, [])) + len(all_times_2.get(4, []))

    num_sample_points = 256
    dt = age_max / num_sample_points

    for age in np.arange(0, age_max + dt, dt):
        toterr = 0
        for vertex_type in [3, 4]:
            s1 = 1
            s2 = 1
            if vertex_type in all_times_1 and len(all_times_1[vertex_type]) > 0:
                times_vector_1 = pyffwd.list_to_vector(all_times_1[vertex_type])
                s1 = pyffwd.surv_kaplan_meier(age, times_vector_1,
                                              reference_pop_1)
            if vertex_type in all_times_2 and len(all_times_2[vertex_type]) > 0:
                times_vector_2 = pyffwd.list_to_vector(all_times_2[vertex_type])
                s2 = pyffwd.surv_kaplan_meier(age, times_vector_2,
                                              reference_pop_2)
            error = (s1 - s2) ** 2
            toterr += error

        toterr = np.sqrt(toterr / 2)
        print(f"{age:.20f},{toterr:.20f},")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Call this program with:")
        print("    python gillespie_errors.py seed runs")
        sys.exit(1)

    seed = int(sys.argv[1])
    runs = int(sys.argv[2])
    main(seed, runs)
