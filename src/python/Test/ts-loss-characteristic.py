import pyffwd
import sys

def main(dt, type):
    # Convert command line arguments
    dt = float(dt)
    type = int(type)

    # System coefficients
    rloh = 5e-7
    mu = 5e-8

    # Create model
    model = pyffwd.Model.create_model(
        5, 
        birth_rates=[0, 0.05, 0.03, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e6, 0, 0, 0, 0]
    )
    model.m_migr = [
        {1: mu, 2: rloh},         # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 4: 0.5 * rloh},  # From vertex 1 to 3 and 4
        {4: 0.5 * mu},            # From vertex 2 to 4
        {},                       # No migration from vertex 3
        {}                        # No migration from vertex 4
    ]

    # Initial q-values
    qvaluesBoth = pyffwd. RealVector([1, 1, 1, 0, 0])
    qvaluesOther = pyffwd. RealVector([1, 1, 1, 0, 0])
    qvaluesOther[type] = 1

    time = 0.0
    tmax = 380.0
    t_write_step = 1.0

    # Print header
    print("age, p1, p2,")

    # Print initial values
    print(f"{time:.17f}, 1.0, 0.0,")
    t_write = time + t_write_step

    while time < tmax:
        while time < t_write - dt:
            pyffwd.heun_q_step(qvaluesBoth, time, dt, model)
            pyffwd.heun_q_step(qvaluesOther, time, dt, model)
            time += dt

        delta = t_write - time
        pyffwd.heun_q_step(qvaluesBoth, time, delta, model)
        pyffwd.heun_q_step(qvaluesOther, time, delta, model)
        time = t_write

        probneither = pyffwd.generating_function(qvaluesBoth,
                                                 model.m_initial_pops)
        probother = pyffwd.generating_function(qvaluesOther,
                                               model.m_initial_pops)
        prob = probneither / probother

        print(f"{time:.17f}, {prob:.17f}, {1.0 - prob:.17f},")
        t_write += t_write_step

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Call this program with\n ./tsconj dt type\n")
        print("type should be 3 or 4")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])
