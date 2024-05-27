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
    plt.savefig(f"P{title[-3:]}.png", transparent=True)

def difference_method(f, a, b, h, w0, w1):
    w_list = [w0, w1]
    t_list = [a, a + h]
    for i in range(int((b - a) / h) - 1):
        w_list.append(4 * w_list[-1] - 3 * w_list[-2] - 2 * h * f(t_list[-2], w_list[-2]))
        t_list.append(a + (i + 2) * h)
    return t_list, w_list

def f(t, y):
    return 1- y

t_values_b, approximation_b = difference_method(f, 0, 1, 0.1, 0, 1 - math.exp(-0.1))
plot_array(t_values_b, approximation_b, "Approximation Using the Difference Method: 16b")

t_values_c, approximation_c = difference_method(f, 0, 1, 0.01, 0, 1 - math.exp(-0.01))
plot_array(t_values_c, approximation_c, "Approximation Using the Difference Method: 16c")