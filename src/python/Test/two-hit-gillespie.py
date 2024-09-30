import threading
import pyffwd

# Number of hardware threads available
num_thr = 1
runs_per_thr = int(1e7)
seed = 1

# System coefficients
mu0 = 0.001
mu1 = 0.001

# Model setup
model = pyffwd.Model.create_model(
        3, 
        birth_rates=[1.0, 1.2, 1.0],
        death_rates=[1.0, 1.0, 1.0],
        initial_pops=[1, 0, 0]
    )
model.m_migr = [
    {1: mu0},
    {2: mu1},
    {}
]

def simulate_times(index):
    # Create a RealVector to hold the results
    results = pyffwd.RealVector()
    
    # Run the simulation
    pyffwd.times_to_final_vertex(model, seed + index, runs_per_thr, 2, results)
    
    # Convert RealVector to Python list
    result_list = list(results)
    return result_list


# Multi-threaded simulation
times = [None] * num_thr
threads = []

for i in range(num_thr):
    thread = threading.Thread(target=lambda idx=i:
                              times.__setitem__(idx, simulate_times(idx)))
    threads.append(thread)
    thread.start()

# Wait for all threads to complete
for thread in threads:
    thread.join()

# Gather the results
all_times = []
for time_list in times:
    if time_list is not None:
        all_times.extend(time_list)

all_times.sort()

# Convert to RealVector
all_times_real_vector = pyffwd.RealVector(all_times)

# Kaplan-Meier plot
pyffwd.print_kaplan_meier(100, all_times_real_vector)

print()
