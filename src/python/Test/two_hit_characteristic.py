import pyffwd

def main():
    # Initialize the model with 3 vertices
    model = pyffwd.Model(3)
    
    # Create RealVector for each attribute
    birth_rates = pyffwd. list_to_vector([1.0, 1.2, 1.0])
    death_rates = pyffwd. list_to_vector([1.0, 1.0, 1.0])
    initial_pops = pyffwd. list_to_vector([1, 0, 0])
    
    # Set the model parameters
    model.m_migr = [{1: 0.001}, {2: 0.001}, {}]
    model.m_birth = birth_rates
    model.m_death = death_rates
    model.m_initial_pops = initial_pops
    
    # We want the probability that site 2 is unoccupied: the corresponding set
    # of q-coordinates is (1,1,0)
    qvalues = pyffwd.RealVector([1, 1, 0])

    time = 0.0
    tmax = 100.0
    dt = 1.0
    print(f"{time}, 1, 0,")

    while time < tmax:
        # Subdivide the time step dt into smaller, finer time steps
        subdivision = 16
        dt2 = dt / subdivision
        for _ in range(subdivision):
            pyffwd.heun_q_step(qvalues, time, dt2, model)

        # Increment time
        time += dt
        
        prob = pyffwd.generating_function(qvalues, model.m_initial_pops)

        print(f"{time}, {prob}, {1.0 - prob},")

if __name__ == "__main__":
    main()
