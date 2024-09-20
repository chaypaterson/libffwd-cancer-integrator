import pyffwd

def main():
    # Create model with 6 stages
    parameters = pyffwd.Model.create_model(
        6, 
        birth_rates= [0, 0, 0.2, 0.27, 0.27, 0],
        death_rates=[0, 0, 0, 0, 0, 0],
        initial_pops=[1e8, 0, 0, 0, 0, 0]
    )
    # System coefficients
    parameters.m_migr = [
        {1: 2.86e-4},
        {2: 1.06e-5},
        {3: 9.00e-7},
        {4: 1.36e-4},
        {5: 4.56e-7}
    ]

    # Initialize qvalues
    default_qvalues = [1, 1, 1, 1, 1, 1]
    qvalues = [pyffwd.RealVector(default_qvalues) 
               for _ in range(parameters.m_stages)]
    
    for vertex in range(parameters.m_stages):
        qvalues[vertex][vertex] = 0

    time = 0.0
    tmax = 80.0
    dt = 1.0

    while time < tmax:
        subdivision = 20
        dt2 = dt / subdivision
        for vertex in range(parameters.m_stages):
            for _ in range(subdivision):
                pyffwd.heun_q_step(qvalues[vertex], time, dt2, parameters)

        time += dt

        # Compute probabilities and print results in the desired format
        output = f"{int(time)}, "  # Formatting the time as an integer
        for vertex in range(parameters.m_stages):
            prob = pyffwd.generating_function(qvalues[vertex],
                                              parameters.m_initial_pops)
            # Limit to 6 significant digits:
            formatted_prob = f"{1.0 - prob:.6g}"
            output += f"{formatted_prob}, "

        print(output)

if __name__ == "__main__":
    main()
