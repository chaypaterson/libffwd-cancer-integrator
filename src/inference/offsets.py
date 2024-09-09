import math

# estimate target -log likelihood values for each study

logl=[10.806779, 59.417208, 238.8904, 902.63374, 3249.0576, 11927.082, 41872.248]
sizes=[10, 32, 100, 320, 1000, 3200, 10000]
offsets=[0.5 * n * math.log(n) for n in sizes]

print(sizes, "\n", offsets, "\n", logl)

print([offsets[n] - logl[n] for n in range(len(sizes))])
print([(offsets[n] - logl[n]) / sizes[n] for n in range(len(sizes))])
print([(offsets[n] - logl[n])**2 / (sizes[n] * sizes[n]) for n in range(len(sizes))])
