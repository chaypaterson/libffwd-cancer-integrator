import pyffwd
import math

EPSILON = 1e-10

def map_onto_data(params, this_data):
    qcorner = pyffwd.RealVector([1] * params.m_stages)
    dt = 0.01
    total = 0

    for age, node in this_data:
        qvals = pyffwd.RealVector(qcorner)
        qvals[node] = 0
        time = 0.0

        while time < age:
            pyffwd.heun_q_step(qvals, time, dt, params)
            time += dt

        prob = pyffwd.generating_function(qvals, params.m_initial_pops)
        pyffwd.heun_q_step(qvals, time, dt, params)
        prob2 = pyffwd.generating_function(qvals, params.m_initial_pops)
        dprob = (prob - prob2) / dt

        dprob_safe = dprob if dprob > EPSILON else EPSILON
        log_dprob = -math.log(dprob_safe)

        if prob < 0.0001:
            prob_str = f"{prob:.7e}"
        else:
            prob_str = f"{prob:.7f}"

        # Format dprob with 7 decimal places to capture full precision
        if dprob_safe < 0.0001:
            dprob_str = f"{dprob_safe:.7e}"
        else:
            dprob_str = f"{dprob_safe:.7f}"

        log_dprob_str = f"{log_dprob:.7f}"

        # Print formatted output
        print(f"{age:.0f}, {node}, {prob_str}, "
              f"{dprob_str}, "
              f"{log_dprob_str}")

        total += log_dprob

    return total

def logsurvival(params, node):
    alpha = params.m_birth[node]
    beta = params.m_death[node]
    if (alpha + beta) != 0:
        psurv = alpha / (alpha + beta)
    else:
        psurv = 1.0
    return -math.log(psurv) if psurv > EPSILON else float('-inf')

def main():
    rloh = 0.5e-2
    mu = 0.5e-3

    params = pyffwd.Model.create_model(
        5, 
        birth_rates= [0, 0.2, 0.2, 0, 0],
        death_rates=[0, 0, 0, 0, 0],
        initial_pops=[1e2, 0, 0, 0, 0]
    )
    params.m_migr = [
        {1: mu, 2: rloh}, # From vertex 0 to 1 and 2
        {3: 0.5 * mu, 4: 0.5 * rloh}, # From vertex 1 to 3 and 4
        {4: 0.5 * mu} # From vertex 2 to 4
    ]

    all_times = [(33, 3), (50, 4)]

    # Perform the calculations and print the results
    total = map_onto_data(params, all_times)
    print(f"{total:.3f}")

if __name__ == "__main__":
    main()
