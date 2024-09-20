import pyffwd

def main():
    # Initialize the model with 3 vertices
    model = pyffwd.Model.create_model(
        3, 
        birth_rates=[1.0, 1.2, 1.0],
        death_rates=[1.0, 1.0, 1.0],
        initial_pops=[1, 0, 0]
    )
    
    # Set the model parameters
    model.m_migr = [{1: 0.001}, {2: 0.001}, {}]
    
    # We want the probability that site 2 is unoccupied: the corresponding set
    # of q-coordinates is (1,1,0)
    qvalues = pyffwd.RealVector([1, 1, 0])

    time = 0.0
    tmax = 100.0
    dt = 1.0
    print(f"{int(time)}, 1, 0,")

    while time < tmax:
        # Subdivide the time step dt into smaller, finer time steps
        subdivision = 16
        dt2 = dt / subdivision
        for _ in range(subdivision):
            pyffwd.heun_q_step(qvalues, time, dt2, model)

        # Increment time
        time += dt
        
        prob = pyffwd.generating_function(qvalues, model.m_initial_pops)

        print(f"{int(time)}, {prob:.6g}, {1.0 - prob:.6g},")

if __name__ == "__main__":
    main()
