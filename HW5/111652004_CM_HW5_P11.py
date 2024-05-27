import math
import matplotlib.pyplot as plt

def plot_arrays(x_axis, y_1s, y_2s, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_2s, 'o', label='Actual Solution')
    plt.plot(x_axis, y_1s, '^', color="black", markerfacecolor='none', label='Approximated Solution by Adams-Moulton Method')

    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-3:]}.png", transparent=True)

def true_y(t):
    return 1 - math.log(1 - math.e * t)

omega_list = [true_y(0.01 * i) for i in range(3)]

def fix_point_iteration(f):
    omega = [omega_list[-j] for j in range(3, 0, -1)]
    new_omega = f(omega_list[-1], omega)
    old_omega = -100
    while abs(old_omega - new_omega) > 1e-8:
        old_omega = new_omega
        new_omega = f(new_omega, omega)
    return new_omega

def g(w, W):
    [w2, w1, w0] = W
    return w0 + 0.01 * (9 * math.exp(w) + 19 * math.exp(w0) - 5 * math.exp(w1) + math.exp(w2)) / 24

for _ in range(3, 21):
    omega_list.append(fix_point_iteration(g))

t_values = [0.01 * i for i in range(21)]
real_y_values = [true_y(0.01 * i) for i in range(21)]
plot_arrays(t_values, omega_list, real_y_values,  "Approximation by Adams-Moulton Method (Functional Iteration): 11a")