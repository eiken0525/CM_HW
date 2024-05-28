import math
import matplotlib.pyplot as plt

def plot_array(x_axis, y_axis, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_axis, 'o-', label='Approximated Solution')
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
    m = len(alpha)

    for i in range(1, N + 1):
        if i < 4:
            k = [[], [], [], []]
            for j in range(m):
                k[0].append(h * f[j](t_values[-1], approx_soln_list[-1]))
            for j in range(m):
                k[1].append(h * f[j](t_values[-1] + h / 2, (approx_soln_list[-1][0] + k[0][0] / 2, approx_soln_list[-1][1] + k[0][1] / 2)))
            for j in range(m):
                k[2].append(h * f[j](t_values[-1] + h / 2, (approx_soln_list[-1][0] + k[1][0] / 2, approx_soln_list[-1][1] + k[1][1] / 2)))
            for j in range(m):
                k[3].append(h * f[j](t_values[-1] + h, (approx_soln_list[-1][0] + k[2][0], approx_soln_list[-1][1] + k[2][1])))
            w = []
            for j in range(m):
                w.append(approx_soln_list[-1][j] + (k[0][j] + 2 * k[1][j] + 2 * k[2][j] + k[3][j]) / 6)
            approx_soln_list.append(w)
            t = a + i * h
            t_values.append(t)
        else:
            w = approx_soln_list[-4:]
            t = t_values[-4:]
            t_values.append(a + i * h)
            w_0_list = []
            W = []
            for j in range(m):
                w_0 = w[3][j] + h * (55 * f[j](t[3], w[3]) - 59 * f[j](t[2], w[2]) + 37 * f[j](t[1], w[1]) - 9 * f[j](t[0], w[0])) / 24
                w_0_list.append(w_0)
            for j in range(m):
                w_0 = w[3][j] + h * (9 * f[j](t_values[-1], w_0_list) + 19 * f[j](t[3], w[3]) - 5 * f[j](t[2], w[2]) + f[j](t[1], w[1])) / 24
                W.append(w_0)
            approx_soln_list.append(W)

    return t_values, approx_soln_list

def phi(t, y):
    return -32.14 / 2 * math.sin(y[0])

def psi(t, y):
    return y[1]

f = [psi, phi]

t_list, y_approx_list = adams_fourth_order_predictor_correction_method(f, (math.pi / 6, 0), 0, 2, 20)

plot_array(t_list, [y[0] for y in y_approx_list], "Approximation by Adams Fourth-Order Predictor-Corrector Algorithm: 15")

for a_y in y_approx_list:
    print(f"{a_y[0]:.6f}")