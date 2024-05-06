import math

TOL = 0.000001
true_a = (math.sin(math.pi ** 2 / 2)) ** 2
true_b = math.pi ** 2 - 4

def cSimpson(f, n, a, b):
    h = (b - a) / n
    def x(i):
        return a + i * h
    acc = f(a) + f(b)
    for i in range(n//2 - 1):
        acc += 2 * f(x(2*(i+1)))
    for i in range(n//2):
        acc += 4 * f(x(2*i+1))
    return h * acc / 3

def adaptiveQuadrature(f, a, b, TOL, num_of_recc):
    def recursive(f, a, b, TOL, fa, fb, fc, approx_0, num_of_recc):
        c = (a + b) / 2; h = (b - a) / 2
        d = (a + c) / 2; e = (c + b) / 2
        fd = f(d); fe = f(e)
        approx_l = h * (fa + 4 * fd + fc) / 6
        approx_r = h * (fc + 4 * fe + fb) / 6
        approx = approx_l + approx_r
        if abs(approx - approx_0) <= TOL:
            num_of_recc[0] += 1
            return approx
        else:
            return (recursive(f, a, c, TOL / 2, fa, fc, fd, approx_l, num_of_recc) +
                    recursive(f, c, b, TOL / 2, fc, fb, fe, approx_r, num_of_recc))

    fa = f(a); fb = f(b); fc = f((a + b) / 2)

    approx_0 = (b - a) * (fa + 4 * fc + fb) / 6

    return (recursive(f, a, b, 10 * TOL, fa, fb, fc, approx_0, num_of_recc), num_of_recc[0])

def f_a(x):
    return x * math.sin(x ** 2)

def f_b(x):
    return x ** 2 * math.sin(x)

num_of_iter_a = 4

print("a.")
former_Simpson_value_a = 0
while(True):
    present_Simpson_value_a = cSimpson(f_a, num_of_iter_a, 0, math.pi)
    if (abs(present_Simpson_value_a - true_a) < TOL and abs(former_Simpson_value_a - true_a) < TOL):
        print(f"It takes {num_of_iter_a} steps to obtain successive approximations with composite Simson's rule for a.")
        break
    else:
        former_Simpson_value_a = present_Simpson_value_a
        num_of_iter_a += 2

num_of_recc_a = [1]
AQ_a, num_of_recc_a = adaptiveQuadrature(f_a, 0, math.pi, TOL, num_of_recc_a)
print(f"For apadative quadrature, it takes {num_of_recc_a} steps.")

num_of_iter_b = 4

print()
print("b.")
former_Simpson_value_b = 0
while(True):
    present_Simpson_value_b = cSimpson(f_b, num_of_iter_b, 0, math.pi)
    if (abs(present_Simpson_value_b - true_b) < TOL and abs(former_Simpson_value_b - true_b) < TOL):
        print(f"It takes {num_of_iter_b} steps to obtain successive approximations with composite Simson's rule for b.")
        break
    else:
        former_Simpson_value_b = present_Simpson_value_b
        num_of_iter_b += 2

num_of_recc_b = [1]
AQ_b, num_of_recc_b = adaptiveQuadrature(f_b, 0, math.pi, TOL, num_of_recc_b)
print(f"For apadative quadrature, it takes {num_of_recc_b} steps.")