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
