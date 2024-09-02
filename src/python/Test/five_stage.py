import pyffwd

def list_to_real_vector(py_list):
    """Helper function to convert a Python list to a std::vector<real_t>."""
    return pyffwd.list_to_vector(py_list)

def main():
    # Create model with 6 stages
    parameters = pyffwd.Model(6)

    # System coefficients
    parameters.m_migr = [
        {1: 2.86e-4},
        {2: 1.06e-5},
        {3: 9.00e-7},
        {4: 1.36e-4},
        {5: 4.56e-7}
    ]
    
    # Convert Python lists to RealVector
    parameters.m_birth = list_to_real_vector([0, 0, 0.2, 0.27, 0.27, 0])
    parameters.m_death = list_to_real_vector([0, 0, 0, 0, 0, 0])
    parameters.m_initial_pops = list_to_real_vector([1e8, 0, 0, 0, 0, 0])

    # Initialize qvalues
    default_qvalues = [1, 1, 1, 1, 1, 1]
    # what is the function of this next line??????? - Chay
    qvalues = [list_to_real_vector(default_qvalues) for _ in range(parameters.m_stages)]
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
        
        # Remove trailing comma and space, then print the formatted output
        print(output.strip(', '))

if __name__ == "__main__":
    main()
