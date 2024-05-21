import math
import matplotlib.pyplot as plt

def plot_arrays(x_axis, y_1s, y_2s, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_1s, 'o', label='Approximated Solution')
    plt.plot(x_axis, y_2s, 'x', label='True Solution')

    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-2:]}.png", transparent=True)

def plot_errors(x_axis, y_1s, y_2s, title):
    plt.figure(figsize=(10, 6))
    y_diff = [abs(y1 - y2) for y1, y2 in zip(y_1s, y_2s)]
    plt.plot(x_axis, y_diff, 's-', label='Difference Between Approximation and True Solution')

    for i in range(len(x_axis)):
        plt.annotate(f'{abs(y_1s[i] - y_2s[i]):.2e}', (x_axis[i], min(y_1s[i], y_2s[i])),
                     textcoords="offset points", xytext=(0,-15), ha='center')

    plt.xlabel('t')
    plt.ylabel('Error')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-2:]}e.png", transparent=True)

def euler_method(f, true_y, alpha, a, b, h):
    N = int((b - a) / h)
    y_0 = alpha
    approx_soln_list = [y_0]
    real_soln_list = [true_y(a)]
    t_values = [a]

    for i in range(1, N + 1):
        t_i = a + i * h
        t_values.append(t_i)
        y_0 += h * f(t_i, y_0)
        approx_soln_list.append(y_0)
        real_soln_list.append(true_y(t_i))

    return t_values, approx_soln_list, real_soln_list

def f_a(t, y):
    return (2 - 2 * t * y) / (t * t + 1)

def true_y_a(t):
    return (2 * t + 1) / (t * t + 1)

def f_b(t, y):
    return y * y / (1 + t)

def true_y_b(t):
    return -1 / math.log(t + 1)

# Get the solutions using Euler's method
t_values_a, approx_soln_a, real_soln_a = euler_method(f_a, true_y_a, 1, 0, 1, 0.1)

# Plot the results
plot_arrays(t_values_a, approx_soln_a, real_soln_a, "Euler's Method Approximation vs. True Solution: 3a")
plot_errors(t_values_a, approx_soln_a, real_soln_a, "Error Plot Graph for 3a")

# Get the solutions using Euler's method
t_values_b, approx_soln_b, real_soln_b = euler_method(f_b, true_y_b, -1 / math.log(2), 1, 2, 0.1)

# Plot the results
plot_arrays(t_values_b, approx_soln_b, real_soln_b, "Euler's Method Approximation vs. True Solution: 3b")
plot_errors(t_values_b, approx_soln_b, real_soln_b, "Error Plot Graph for 3b")