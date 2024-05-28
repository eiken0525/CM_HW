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

def taylor_method(n, f_family, real_y, alpha, a, b, h):
    N = int((b - a) / h)
    y_0 = alpha
    approx_soln_list = [y_0]
    real_soln_list = [y_0]
    t_values = [a]

    for i in range(1, N + 1):
        T = 0
        for ii in range(n):
            T += h ** ii * f_family[ii](t_values[-1], approx_soln_list[-1]) / math.factorial(ii + 1)
        approx_soln_list.append(approx_soln_list[-1] + h * T)
        t_values.append(a + h * i)
        real_soln_list.append(real_y(a + h * i))

    return t_values, approx_soln_list, real_soln_list

def f(t, y):
    return 2 * y / t + t ** 2 * math.exp(t)

def Df(t, y):
    return 2 * y / t ** 2 + (4 * t + t ** 2) * math.exp(t)

def D2f(t, y):
    return (t ** 2 + 6 * t + 6) * math.exp(t)

def D3f(t, y):
    return (t ** 2 + 8 * t + 12) * math.exp(t)

def y(t):
    return t ** 2 * (math.exp(t) - math.e)

f_family = [f, Df, D2f, D3f]


t_values_a, approx_soln_a, real_soln_a = taylor_method(2, f_family, y, 0, 1, 2, 0.1)

plot_arrays(t_values_a, approx_soln_a, real_soln_a, "Taylor's Method of Order Two Approximated Solution: 6a")

for i, j, k in zip(t_values_a, approx_soln_a, real_soln_a):
    print(f"${i :.1f}$ & ${j :.10f}$ & ${k :.10f}$ \\\\")
    print("\hline")

t_values_c, approx_soln_c, real_soln_c = taylor_method(4, f_family, y, 0, 1, 2, 0.1)

plot_arrays(t_values_c, approx_soln_c, real_soln_c, "Taylor's Method of Order Four Approximated Solution: 6c")

for i, j, k in zip(t_values_c, approx_soln_c, real_soln_c):
    print(f"${i :.1f}$ & ${j :.10f}$ & ${k :.10f}$ \\\\")
    print("\hline")

for soln_list in [approx_soln_a]:
    print(f"{0.6 * soln_list[0] + 0.4 * soln_list[1] :.10f}")
    print(f"{0.5 * soln_list[5] + 0.5 * soln_list[6] :.10f}")
    print(f"{0.3 * soln_list[9] + 0.7 * soln_list[10] :.10f}")

for t in [1.04, 1.55, 1.97]:
    print(f"{y(t):.10f}")

def cubic_hermite_interpolation(tuple1, tuple2, x):
        x1, fx1, dfx1 = tuple1
        x2, fx2, dfx2 = tuple2

        df = (fx2-fx1)/(x2-x1)
        ddf1 = (df - dfx1)/(x2-x1)
        ddf2 = (dfx2 - df)/(x2-x1)
        dddf = (ddf2 - ddf1)/(x2-x1)

        return fx1 + (x-x1)*dfx1 + (x-x1)**2 * ddf1 + (x-x1)**2 * (x-x2)*dddf

print(f"1.04: {cubic_hermite_interpolation((1, 0, f(1, 0)), (1.1, approx_soln_c[1], f(1.1, approx_soln_c[1])), 1.04): .10f}")
print(f"1.55: {cubic_hermite_interpolation((1.5, approx_soln_c[5], f(1.5, approx_soln_c[5])), \
    (1.6, approx_soln_c[6], f(1.6, approx_soln_c[6])), 1.55): .10f}")
print(f"1.97: {cubic_hermite_interpolation((1.9, approx_soln_c[9], f(1.9, approx_soln_c[9])), \
    (2, approx_soln_c[10], f(2, approx_soln_c[10])), 1.97): .10f}")