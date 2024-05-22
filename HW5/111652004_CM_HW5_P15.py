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

def adams_fourth_order_predictor_correction_method(f, alpha, a, b, N):
    h = (b - a) / N
    t_0 = a
    y_0 = alpha
    approx_soln_list = [y_0]
    t_values = [t_0]

    for i in range(1, N + 1):
        if i < 4:
            t = t_values[-4]
            w = approx_soln_list[-4]
            k_1 = h * f(t[-1], w)
            k_2 = h * f(t[-1] + h / 2, w + k_1 / 2)
            k_3 = h * f(t[-1] + h / 2, w + k_2 / 2)
            k_4 = h * f(t[-1] + h, w + k_3)
            approx_soln_list.append(w + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6)
            t_values.append(a + i * h)
        else:
            t_values.append(a + i * h)
            w_0 = w[0] + h * (55 * f(t[0], w[0]) - 59 * f(t[1], w[1]) + 37 * f(t[2], w[2]) - 9 * f(t[3], w[3])) / 24
            w_0 = w[0] + h * (9 * f(t_values[-1], w_0) + 19 * f(t[0], w[0]) - 5 * f(t[1], w[1]) + f(t[2], w[2])) / 24
            approx_soln_list.append(w_0)

    return t_values, approx_soln_list

def f(t, y):
    return 0
