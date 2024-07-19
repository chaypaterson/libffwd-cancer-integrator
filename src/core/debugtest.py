import pybinding


# Create an instance of the Model class
model = pybinding.Model(3)
model.m_stages = 3
model.m_birth = pybinding.list_to_vector([0.1, 0.2, 0.3])
model.m_death = pybinding.list_to_vector([0.05, 0.05, 0.05])
#model.m_birth([0.1, 0.2, 0.3]) TODO
#model.m_death([0.05, 0.05, 0.05]) TODO
model.m_migr = [
    {1: 0.01, 2: 0.02},  # From vertex 0 to vertices 1 and 2
    {0: 0.01, 2: 0.01},  # From vertex 1 to vertices 0 and 2
    {0: 0.02, 1: 0.01}   # From vertex 2 to vertices 0 and 1
]

# TODO 4. check members are mutable (lower priority)

# Define initial q-coordinates and time step
qcoords = [0.5, 0.5, 0.5]

# Convert list to std::vector<real_t>
vec_qcoords = pybinding.list_to_vector(qcoords)

# TODO HIGH PRIORITY:
# TODO export vector<real_t> instead of converting from lists
# TODO: 1. bindings for vectors
# TODO: 2. initialisers for vectors, converter from lists/dictionaries to vectors.

time = 0.0
dt = 1.0

# TODO 3. struct/class containing these three? (lower priority)

# Test zeroing function:
pybinding.test_reference(vec_qcoords)
print("Values after test_reference:", vec_qcoords)

# Test rhs_flow function
result_rhs_flow = pybinding.rhs_flow(vec_qcoords, model)
print("Result from rhs_flow:", result_rhs_flow)

# Perform Heun's method step
qcoords = [0.5, 0.5, 0.5]
vec_qcoords = pybinding.list_to_vector(qcoords)
pybinding.heun_q_step(vec_qcoords, time, dt, model)
print("Values after one Heun's method step:", vec_qcoords)

# Perform Implicit Euler method step
qcoords = [0.5, 0.5, 0.5]
vec_qcoords = pybinding.list_to_vector(qcoords)
pybinding.implicit_q_step(vec_qcoords, time, dt, model)
print("Values after Implicit Euler method step:", vec_qcoords)
