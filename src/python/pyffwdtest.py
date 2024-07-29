import pyffwd
import gsl
import random

# Create an instance of the Model class
model = pyffwd.Model(3)
model.m_stages = 3
model.m_birth = pyffwd.list_to_vector([0.1, 0.2, 0.3])
model.m_death = pyffwd.list_to_vector([0.05, 0.05, 0.05])
#model.m_birth([0.1, 0.2, 0.3])
#model.m_death([0.05, 0.05, 0.05])
model.m_migr = [
    {1: 0.01, 2: 0.02},  # From vertex 0 to vertices 1 and 2
    {0: 0.01, 2: 0.01},  # From vertex 1 to vertices 0 and 2
    {0: 0.02, 1: 0.01}   # From vertex 2 to vertices 0 and 1
]

# Define initial q-coordinates and time step
qcoords = [0.5, 0.5, 0.5]
vec_qcoords = pyffwd.list_to_vector(qcoords)
time = 0.0
dt = 1.0

# Perform rhs_flow function
result_rhs_flow = pyffwd.rhs_flow(vec_qcoords, model)
print("Result from rhs_flow:", result_rhs_flow)

# Perform Heun's method step
qcoords = [0.5, 0.5, 0.5]
vec_qcoords = pyffwd.list_to_vector(qcoords)
pyffwd.heun_q_step(vec_qcoords, time, dt, model)
print("Values after one Heun's method step:", vec_qcoords)

# Perform Implicit Euler method step
vec_qcoords = pyffwd.list_to_vector(qcoords)
pyffwd.implicit_q_step(vec_qcoords, time, dt, model)
print("Values after Implicit Euler method step:", vec_qcoords)

# Perform rungekutta method
vec_qcoords = pyffwd.list_to_vector(qcoords)
pyffwd.rungekutta_q_step(vec_qcoords, time, dt, model)
print("Values after rungekutta method step:", vec_qcoords)

# Calculate generating function
initial_pops = [10, 20, 30]
vec_qcoords = pyffwd.list_to_vector(qcoords)
initial_pops = pyffwd.list_to_vector(initial_pops)
psi = pyffwd.generating_function(vec_qcoords, initial_pops)
print("Generating Function:", psi)

# Perform first_passage_time_single & first_passage_time_multiple
# Create an RNG instance
rng = pyffwd.RNGWrapper(seed=42)
final_vertex = 2
# Call the first_passage_time_single method on the RNGWrapper instance
time_to_final_vertex = rng.first_passage_time_single(model, final_vertex)
print(f"First passage time to final vertex {final_vertex}: {time_to_final_vertex}")
# Define multiple final vertices for which you want to calculate the first passage time
final_vertices = [1, 2]
results = rng.first_passage_time_multiple(model, final_vertices)
print(f"First passage times to final vertices {final_vertices}: {results}")

# Perform print_results
print("print_results...")
all_times = [random.uniform(0, 10) for _ in range(10)]
all_times_vector = pyffwd.list_to_vector(all_times)
pyffwd.print_results(all_times_vector)

# Perform print_kaplan_meier
print("print_kaplan_meier...")
pyffwd.print_kaplan_meier(10.0, all_times_vector)
print("print_kaplan_meier with ref_pop...")
pyffwd.print_kaplan_meier(10.0, all_times_vector, ref_pop=5)

# Perform surv_kaplan_meier
print("surv_kaplan_meier:")
survival_rate = pyffwd.surv_kaplan_meier(5.0, all_times_vector, ref_pop=5)
print(survival_rate)
assert isinstance(survival_rate, float)

# Perform print_naive_estimator
print("print_naive_estimator...")
pyffwd.print_naive_estimator(10.0, all_times_vector)
