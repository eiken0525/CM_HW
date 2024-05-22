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

def runge_kutta_fehlberg_method(f, alpha, TOL, a, b, hmax, hmin):
    h = hmax
    y_0 = alpha
    approx_soln_list = [y_0]
    t_values = [a]

    FLAG = True

    while(FLAG):
        t = t_values[-1]
        k_1 = h * f(t, y_0)
        k_2 = h * f(t + h / 4, y_0 + k_1 / 4)
        k_3 = h * f(t + 3 * h / 8, y_0 + 3 * k_1 / 32 + 9 * k_2 / 32)
        k_4 = h * f(t + 12 * h / 13, y_0 + 1932 * k_1 / 2197 - 7200 * k_2 / 2197 + 7296 * k_3 / 2197)
        k_5 = h * f(t + h, y_0 + 439 * k_1 / 216 - 8 * k_2 + 3680 * k_3 / 513 - 845 * k_4 / 4104)
        k_6 = h * f(t + h / 2, y_0 - 8 * k_1 / 27 + 2 * k_2 - 3544 * k_3 / 2565 + 1859 * k_4 / 4104 - 11 * k_5 / 40)

        R = abs(k_1 / 360 - 128 * k_3 / 4275 - 2197 * k_4 / 75240 + k_5 / 50 + 2 * k_6 / 55) / h

        if R <= TOL:
            t += h
            y_0 += 25 * k_1 / 216 + 1408 * k_3 / 2565 + 2197 * k_4 / 4104 - k_5 / 5.
            t_values.append(t)
            approx_soln_list.append(y_0)

        delta = 0.84 * pow(TOL / R, 0.25)
        if delta <= 0.1:
            h = 0.1 * h
        elif delta >= 4:
            h = 4 * h
        else:
            h = 8 * h

        if h > hmax:
            h = hmax

        if t >= b:
            FLAG = False
        elif t + h > b:
            h = b - t
        elif h < hmin:
            FLAG = False
            print("Minimum h exceeded.")

    return t_values, approx_soln_list

def f_a(t, y):
    return y * y / (t * t) + t / t

def f_b(t, y):
    return math.sin(t) + math.exp(-t)

t_values_a, approx_soln_a = runge_kutta_fehlberg_method(f_a, 1, 10 ** -4, 1, 1.2, 0.05, 0.02)

plot_array(t_values_a, approx_soln_a, "Runge-Kutta-Fehlberg Method Approximated Solution: 9a")

t_values_b, approx_soln_b = runge_kutta_fehlberg_method(f_b, 0, 10 ** -4, 0, 1, 0.25, 0.02)

plot_array(t_values_b, approx_soln_b, "Runge-Kutta-Fehlberg Method Approximated Solution: 9b")
