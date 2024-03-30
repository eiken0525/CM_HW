def f(i):
    return 135000*i-1000*(1-(1+i)**(-360))

# Set endpoints
a = 0.001
b = 0.01

# Set tolerance to 10^-9
TOL = 0.000000001

# Set the maximum number of iterations
N_0 = 500

FA = f(a)
for n in range(N_0):
    i = a + (b - a) / 2
    FI = f(i)
    if FI == 0 or (b - a) / 2 < TOL:
        print(f"The solution is approximately i_{n + 1} = {i}.")
        exit(1)
    if FA * FI > 0:
        a = i
        FA = FI
    else:
        b = i
        FB = FI
print(f"The bisection method failed after {N_0} times.")


