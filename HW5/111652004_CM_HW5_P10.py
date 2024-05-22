import math
import matplotlib.pyplot as plt

def plot_arrays(x_axis, y_1s, y_2s, y_3s, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_1s, 'o', label='Approximated Solution by Adams-Bashforth Method')
    plt.plot(x_axis, y_2s, 'D', label='Approximated Solution by Runge-Kutta Method')
    plt.plot(x_axis, y_3s, 'x', label='Actual Solution')

    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-3:]}.png", transparent=True)

def adams_bashforth_method(f, true_y, alpha, a, b, h):
    N = int((b - a) / h)
    approx_soln_list = alpha.copy()
    real_soln_list = alpha.copy()
    t_values = [a, a + h, a + 2 * h, a + 3 * h]

    for i in range(4, N + 1):
        w = approx_soln_list[-4:]
        t = t_values[-4:]
        approx_soln_list.append(w[3] + h * (55 * f(t[3], w[3]) - 59 * f(t[2], w[2]) + 37 * f(t[1], w[1]) - 9 * f(t[0], w[0])) / 24)
        t_i = a + i * h
        t_values.append(t_i)
        real_soln_list.append(true_y(t_i))

    return t_values, approx_soln_list, real_soln_list

def runge_kutta_method(f, true_y, alpha, a, b, h):
    N = int((b - a) / h)
    approx_soln_list = [alpha]
    real_soln_list = [alpha]
    t_values = [a]

    for i in range(1, N + 1):
        t_i = a + i * h
        real_soln_list.append(true_y(t_i))
        t_values.append(t_i)
        if i < 1:
            approx_soln_list.append(true_y(t_i))
        else:
            t = t_values[-1]
            w = approx_soln_list[-1]
            k_1 = h * f(t, w)
            k_2 = h * f(t + h / 2, w + k_1 / 2)
            k_3 = h * f(t + h / 2, w + k_2 / 2)
            k_4 = h * f(t_i, w + k_3)
            approx_soln_list.append(w + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6)

    return t_values, approx_soln_list, real_soln_list

def f_a(t, y):
    return y / t - y * y / (t * t)

def true_y_a(t):
    return t / (1 + math.log(t))

def f_b(t, y):
    return -5 * y + 5 * t * t + 2 * t

def true_y_b(t):
    return t * t + math.exp(-5 * t) / 3

t_values_a, approx_soln_a_RK, real_soln_a = runge_kutta_method(f_a, true_y_a, 1, 1, 2, 0.1)

t_values_a, approx_soln_a_AB, real_soln_a = adams_bashforth_method(f_a, true_y_a, approx_soln_a_RK[0:4], 1, 2, 0.1)

plot_arrays(t_values_a, approx_soln_a_AB, approx_soln_a_RK, real_soln_a, "Adams-Bashforth and Runge-Kutta Approximations vs. True Solution: 10a")

t_values_b, approx_soln_b_RK, real_soln_b = runge_kutta_method(f_b, true_y_b, 1 / 3, 0, 1, 0.1)

t_values_b, approx_soln_b_AB, real_soln_b = adams_bashforth_method(f_b, true_y_b, approx_soln_b_RK[0:4], 0, 1, 0.1)

plot_arrays(t_values_b, approx_soln_b_AB, approx_soln_b_RK, real_soln_b, "Adams-Bashforth and Runge-Kutta Approximations vs. True Solution: 10b")