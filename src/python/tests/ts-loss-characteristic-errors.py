import sys
import numpy as np
import pyffwd

def main(time_step):
    # System coefficients
    rloh = 0.5e-2
    mu = 0.5e-3

    # Initialize the Model
    model = pyffwd.Model.create_model(
        5, 
        birth_rates=[0, 0.2, 0.2, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e2, 0, 0, 0, 0]
    )
    
    model.m_migr = [
        {1: mu, 2: rloh},
        {3: 0.5 * mu, 4: 0.5 * rloh},
        {4: 0.5 * mu},
        {},
        {}
    ]

    qvaluesF = pyffwd.RealVector([1, 1, 1, 0, 0])
    qvaluesF2 = pyffwd.RealVector([1, 1, 1, 0, 0])
    qvalues4 = pyffwd.RealVector([1, 1, 1, 0, 1])
    qvalues42 = pyffwd.RealVector([1, 1, 1, 0, 1])

    time = 0.0
    tmax = 100.0
    dt = float(time_step)  # integration step
    half = 0.5 * dt
    t_write_step = 1.0  # write out step

    # Print header
    print("age, err,")
    print(f"{time:.20f}, 0.00000000000000000000, ")

    t_write = time + t_write_step

    while time < tmax:
        while time < t_write - dt:
            pyffwd.heun_q_step(qvaluesF, time, dt, model)
            pyffwd.heun_q_step(qvalues4, time, dt, model)
            # Advance qvalues.2 by half the relevant time step, twice
            for _ in range(2):
                pyffwd.heun_q_step(qvaluesF2, time, half, model)
                pyffwd.heun_q_step(qvalues42, time, half, model)
            time += dt

        delta = t_write - time
        pyffwd.heun_q_step(qvaluesF, time, delta, model)
        pyffwd.heun_q_step(qvalues4, time, delta, model)
        halfdelta = 0.5 * delta
        for _ in range(2):
            pyffwd.heun_q_step(qvaluesF2, time, halfdelta, model)
            pyffwd.heun_q_step(qvalues42, time, halfdelta, model)

        time = t_write

        # Compute the P(t) values
        psiF = pyffwd.generating_function(qvaluesF, model.m_initial_pops)
        psi4 = pyffwd.generating_function(qvalues4, model.m_initial_pops)
        prob = psiF / psi4
        psiF2 = pyffwd.generating_function(qvaluesF2, model.m_initial_pops)
        psi42 = pyffwd.generating_function(qvalues42, model.m_initial_pops)
        prob2 = psiF2 / psi42

        # Compute the corresponding error
        err = (prob - prob2)
        err /= (2 << 2) - 1  # Richardson extrapolation
        err *= (1 - 2 * (err < 0))  # abs value

        # Write out errors
        print(f"{time:.20f}, {err:.20f}")

        t_write += t_write_step

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("please provide a time step: ")
        print("  python this_program.py time_step")
        sys.exit(1)

    time_step = sys.argv[1]
    main(time_step)
