import math
import matplotlib.pyplot as plt

def plot_array(x_axis, y_axis, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_axis, 'o', label='Approximated Solution')
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-2:]}.png", transparent=True)

def taylor_method(n, f_family, alpha, a, b, h):
    N = int((b - a) / h)
    y_0 = alpha
    approx_soln_list = [y_0]
    t_values = [a]

    for i in range(1, N + 1):
        T = 0
        for ii in range(n):
            T += h ** ii * f_family[ii](t_values[-1], approx_soln_list[-1]) / math.factorial(ii + 1)
        approx_soln_list.append(approx_soln_list[-1] + h * T)
        t_values.append(a + h * i)

    return t_values, approx_soln_list

def f_a(t, y):
    return (1 + t) / (1 + y)

def f_b(t, y):
    return -y + t * y ** 0.5

def Df_a(t, y):
    return 1 / (1 + y) - (1 + t) ** 2 / (1 + y) ** 3

def Df_b(t, y):
    return y + y * 0.5 - 3 * t * y ** 0.5 / 2 + t ** 2 / 2

f_a_family = [f_a, Df_a]
f_b_family = [f_b, Df_b]


t_values_a, approx_soln_a = taylor_method(2, f_a_family, 2, 1, 2, 0.5)

plot_array(t_values_a, approx_soln_a, "Taylor's Method of Order Two Approximated Solution: 5a")

t_values_b, approx_soln_b = taylor_method(2, f_b_family, 2, 2, 3, 0.25)

plot_array(t_values_b, approx_soln_b, "Taylor's Method of Order Two Approximated Solution: 5b")