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

def adams_variable_step_size_predictor_corrector(f, alpha, a, b, TOL, hmax, hmin):
    def RK4(h, v0, x0):
        x_v_tuple_list = [(x0, v0)]
        for j in range(1, 4):
            x, v = x_v_tuple_list[-1]
            k1 = h * f(x, v)
            k2 = h * f(x + h / 2, v + k1 / 2)
            k3 = h * f(x + h / 2, v + k2 / 2)
            k4 = h * f(x + h, v + k3)
            V = v + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            X = x0 + j * h
            x_v_tuple_list.append((X, V))
        return x_v_tuple_list
    
    t0, w0, h, FLAG, LAST = a, alpha, hmax, True, False
    t_values = [t0]
    approx_soln_list = [w0]
    RK4_approx = RK4(h, w0, t0)
    NFLAG = True
    i = 4
    t = RK4_approx[-1][0] + h

    while FLAG:
        w_list = [RK4_approx[-j][1] for j in range(4, 0, -1)]
        t_list = [RK4_approx[-j][0] for j in range(4, 0, -1)]

        WP = w_list[-1] + h * (55 * f(t_list[-1], w_list[-1]) - 59 * f(t_list[-2], w_list[-2]) + 37 * f(t_list[-3], w_list[-3]) - 9 * f(t_list[-4], w_list[-4])) / 24
        WC = w_list[-1] + h * (9 * f(t, WP) + 19 * f(t_list[-1], w_list[-1]) - 5 * f(t_list[-2], w_list[-2]) + f(t_list[-3], w_list[-3])) / 24
        sigma = 19 * abs(WC - WP) / (270 * h)
        
        if sigma <= TOL:
            w_list.append(WC)
            t_list.append(t)

            if NFLAG:
                for j in range(i - 3, i + 1):
                    t_values.append(t_list[j - (i - 4)])
                    approx_soln_list.append(w_list[j - (i - 4)])
            else:
                t_values.append(t)
                approx_soln_list.append(WC)

            if LAST:
                FLAG = False
            else:
                i += 1
                NFLAG = False
                if sigma <= 0.1 * TOL or t + h > b:
                    q = (TOL / 2 / sigma) ** 0.25
                    if q > 4:
                        h *= 4
                    else:
                        h *= q
                    if h > hmax:
                        h = hmax
                    if t + 4 * h > b:
                        h = (b - t) / 4
                        LAST = True
                    RK4_approx = RK4(h, w_list[-1], t_list[-1])
                    NFLAG = True
                    i += 3
        else:
            q = (TOL / (2 * sigma)) ** 0.25
            if q < 0.1:
                h *= 0.1
            else:
                h *= q
            if h < hmin:
                FLAG = 0
                raise ValueError("Minimum step size exceeded")
            else:
                if NFLAG:
                    i -= 3
                    RK4_approx = RK4(h, w_list[-1], t_list[-1])
                    NFLAG = True

        t = t_list[-1] + h

    return t_values, approx_soln_list


def f_a(t, y):
    return math.sin(t) + math.exp(-t)

def f_b(t, y):
    return - t * y + 4 * t / y

t_values_a, approximation_a = adams_variable_step_size_predictor_corrector(f_a, 0, 0, 1, 1e-4, 0.2, 0.01)

plot_array(t_values_a, approximation_a, "Approximation by the Adams Variable Step-Siza Predictor-Corrector Algorithm: 13a")

t_values_b, approximation_b = adams_variable_step_size_predictor_corrector(f_b, 1, 0, 1, 1e-4, 0.2, 0.01)

plot_array(t_values_b, approximation_b, "Approximation by the Adams Variable Step-Siza Predictor-Corrector Algorithm: 13b")