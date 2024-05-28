import math
import matplotlib.pyplot as plt

def plot_array(x_axis, y_axis, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_axis, '.', label='$y=P(x)$')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-3:]}.png", transparent=True)

def f(x):
    return x ** 3 + 1.5 * x ** 2 - 3 * x - 0.5

x_values = [-3 + 0.01 * i for i in range(601)]
y_values = [f(-3 + 0.01 * i) for i in range(601)]

plot_array(x_values, y_values, "The graph of $y=\lambda^3+\dfrac{3}{2}\lambda^2-3\lambda-\dfrac{1}{2}$: 17b")