def f(x):
    return x ** 3 - 5 * x ** 2 + 8 * x - 6

# Set initial approximation
x_0 = 4

# Set tolerance
TOL = 0.00001

# Set the maximum number of iterations
N_0 = 500

for i in range(N_0):
    x = x_0 - f(x_0) / 16
    if abs(x - x_0) < TOL:
        print(f"The root of P(x) calculated by Newton's method is approixmately {x}.")
        exit(0)
    x_0 = x

print(f"Newton's method failed after {N_0} iterations.")



