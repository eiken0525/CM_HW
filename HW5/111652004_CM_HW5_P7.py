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

    plt.xlabel('t')
    plt.ylabel('Error')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-2:]}e.png", transparent=True)

def modified_euler_method(f, true_y, alpha, a, b, h):
    N = int((b - a) / h)
    y_0 = alpha
    approx_soln_list = [y_0]
    real_soln_list = [true_y(a)]
    t_values = [a]

    for i in range(1, N + 1):
        t_i = a + i * h
        t_values.append(t_i)
        y_0 += h * (f(t_values[-2], y_0) + f(t_i, y_0 + h * f(t_values[-2], y_0))) / 2
        approx_soln_list.append(y_0)
        real_soln_list.append(true_y(t_i))

    return t_values, approx_soln_list, real_soln_list

def f_a(t, y):
    return y * y / (1 + t)

def true_y_a(t):
    return -1 / math.log(t + 1)

def f_b(t, y):
    return (y * y + y) / t

def true_y_b(t):
    return 2 * t / (1 - 2 * t)


t_values_a, approx_soln_a, real_soln_a = modified_euler_method(f_a, true_y_a, -1 / math.log(2), 1, 2, 0.1)

plot_arrays(t_values_a, approx_soln_a, real_soln_a, "Modified Euler Method Approximation vs. True Solution: 7a")
plot_errors(t_values_a, approx_soln_a, real_soln_a, "Error Plot Graph for 7a")

t_values_b, approx_soln_b, real_soln_b = modified_euler_method(f_b, true_y_b, -2, 1, 3, 0.2)

plot_arrays(t_values_b, approx_soln_b, real_soln_b, "Modified Euler Method Approximation vs. True Solution: 7b")
plot_errors(t_values_b, approx_soln_b, real_soln_b, "Error Plot Graph for 7b")