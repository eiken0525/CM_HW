def P(x):
    return x ** 3 - 5 * x ** 2 + 8 * x - 6

# Set approximations
p_0 = 0; p_1 = 1; p_2 = 2

# Set tolerance
TOL = 0.00001

# Set the maximum number of iterations
N_0 = 500

h_1 = p_1 - p_0; h_2 = p_2 - p_1
delta_1 = (P(p_1) - P(p_0)) / h_1; delta_2 = (P(p_2) - P(p_1)) / h_2
d = (delta_2 - delta_1) / (h_2 + h_1)

for i in range(N_0):
    b = delta_2 + h_2 * d
    D = (b ** 2 - 4 * P(p_2) * d) ** 0.5
    if abs(b - D) < abs(b + D):
        E = b + D
    else:
        E = b - D
    h = -2 * P(p_2) / E
    p = p_2 + h
    if abs(h) < TOL:
        print(f"The solution calculated by Muller's method is approximately {p}.")
        exit(0)
    p_0 = p_1; p_1 = p_2; p_2 = p
    h_1 = p_1 - p_0; h_2 = p_2 - p_1
    delta_1 = (P(p_1) - P(p_0)) / h_1; delta_2 = (P(p_2) - P(p_1)) / h_2
    d = (delta_2 - delta_1) / (h_2 + h_1)

print(f"Muller's method failed after {N_0} iterations.")



