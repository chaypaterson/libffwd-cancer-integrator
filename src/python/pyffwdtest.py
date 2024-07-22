import pyffwd


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