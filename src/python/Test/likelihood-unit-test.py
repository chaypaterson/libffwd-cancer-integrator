import pyffwd
import math

EPSILON = 1e-10

def map_onto_data(params, this_data):
    qcorner = pyffwd.list_to_vector([1] * params.m_stages)
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

        # Output the result in the required format
        dprob_safe = dprob if dprob > EPSILON else EPSILON
        log_dprob = -math.log(dprob_safe)
        print(f"{age}, {node}, {prob:.6f}, {dprob:.6f}, {log_dprob:.5f}")

        total += -log_dprob

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

    params = pyffwd.Model(5)
    params.m_migr[0][1] = mu
    params.m_migr[0][2] = rloh
    params.m_migr[1][3] = 0.5 * mu
    params.m_migr[1][4] = 0.5 * rloh
    params.m_migr[2][4] = 0.5 * mu

    params.m_birth = pyffwd.list_to_vector([0.1, 0.2, 0.2, 0.1, 0.1])
    params.m_death = pyffwd.list_to_vector([0.1, 0.1, 0.1, 0.1, 0.1])
    params.m_initial_pops = pyffwd.list_to_vector([1e2, 1e2, 1e2, 1e2, 1e2])

    all_times = [(33.0, 3), (50.0, 4)]

    # Perform the calculations and print the results
    total = map_onto_data(params, all_times)
    print(f"{total:.3f}")

if __name__ == "__main__":
    main()
