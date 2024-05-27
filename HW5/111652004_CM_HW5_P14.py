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

def extrapolation_method(f, alpha, a, b, TOL, hmax, hmin):
    NK = [0, 2, 4, 6, 8, 12, 16, 24, 32]
    t_values = []
    approx_soln_list = []
    Q = [[0] for _ in range(8)]
    TO, WO, h, FLAG = a, alpha, hmax, True
    for i in range(1, 8):
        for j in range(1, i + 1):
            Q[i].append((NK[i + 1] / NK[j]) ** 2)
    while (FLAG):
        y = [0]
        k, NFLAG = 1, False
        while(k <= 8 and (not NFLAG)):
            HK = h / NK[k]
            T = TO
            W2 = WO
            W3 = W2 + HK * f(T, W2)
            T = TO + HK
            for j in range(1, NK[k]):
                W1 = W2
                W2 = W3
                W3 = W1 + 2 * HK * f(T, W2)
                T = TO + (j + 1) * HK

            y.append((W3 + W2 + HK * f(T, W3)) / 2)
            if k >= 2:
                j = k
                v = y[1]
                while(j >= 2):
                    y[j - 1] = y[j] + (y[j] - y[j - 1]) / (Q[k - 1][j - 1] - 1)
                    j -= 1
                if abs(y[1] - v) <= TOL:
                    NFLAG = True
            k += 1
        k -= 1
        if not NFLAG:
            h /= 2
            if h < hmin:
                exit(1)
        else:
            WO, TO = y[1], TO + h
            t_values.append(TO)
            approx_soln_list.append(WO)
            if TO >= b:
                FLAG = False
            elif TO + h > b:
                h = b - TO
            elif k <= 3 and h < 0.5 * hmax:
                h *= 2
    return t_values, approx_soln_list


def f(t, y):
    return 2.9e-2 * y - 1.4e-7 * y ** 2

t_values, approx_soln = extrapolation_method(f, 50976, 0, 5, 1e-9, 0.1000000001, 0.999999999)

plot_array(t_values, approx_soln, "Approximation to the Logistic Equation: 14")