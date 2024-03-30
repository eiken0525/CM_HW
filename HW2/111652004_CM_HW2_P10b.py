import math

def g(p):
    return math.sqrt(1 + 1 / p)

# Set initial approximation
p_0 = 1.3

# Set tolerance to 10^-4
TOL = 0.0001

# Set the maximum number of iteration
N_0 = 500

for i in range(N_0):
    p_1 = g(p_0)
    p_2 = g(p_1)
    p = p_0 - (p_1 - p_0) ** 2 / (p_2 - 2 * p_1 + p_0)
    if abs(p - p_0) < TOL:
        print(f"The solution calculated by Steffensen's method is approximately {p}.")
        exit(0)
    p_0 = p
print(f"Steffensen's method failed after {N_0} times.")


