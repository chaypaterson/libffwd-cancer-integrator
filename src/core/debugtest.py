import pybinding

# Create an instance of the Model class
model = pybinding.Model(3)
model.m_stages = 3
model.m_birth = [0.1, 0.2, 0.3]
model.m_death = [0.05, 0.05, 0.05]
model.m_migr = [
    {1: 0.01, 2: 0.02},  # From vertex 0 to vertices 1 and 2
    {0: 0.01, 2: 0.01},  # From vertex 1 to vertices 0 and 2
    {0: 0.02, 1: 0.01}   # From vertex 2 to vertices 0 and 1
]

# Define initial q-coordinates and time step
qcoords = [0.5, 0.5, 0.5]
time = 0.0
dt = 1.0

# Test rhs_flow function
result_rhs_flow = pybinding.rhs_flow(qcoords, model)
print("Result from rhs_flow:", result_rhs_flow)

print("Initial values:", qcoords)

# Perform Heun's method step
pybinding.heun_q_step(qcoords, time, dt, model)
print("Values after one Heun's method step:", qcoords)

# Reset qcoords for implicit Euler method
qcoords = [0.5, 0.5, 0.5]
print("Initial values reset:", qcoords)

# Perform Implicit Euler method step
pybinding.implicit_q_step(qcoords, time, dt, model)
print("Values after Implicit Euler method step:", qcoords)

pybinding.test_reference(qcoords)
print("Values after test_reference:", qcoords)



