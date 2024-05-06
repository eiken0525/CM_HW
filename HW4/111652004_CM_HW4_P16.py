import math

TOL = 1e-6

def adaptiveQuadrature(f, a, b, TOL):
    def recursive(f, a, b, TOL, fa, fb, fc, approx_0):
        c = (a + b) / 2; h = (b - a) / 2
        d = (a + c) / 2; e = (c + b) / 2
        fd = f(d); fe = f(e)
        approx_l = h * (fa + 4 * fd + fc) / 6
        approx_r = h * (fc + 4 * fe + fb) / 6
        approx = approx_l + approx_r
        if abs(approx - approx_0) <= TOL:
            return approx
        else:
            return (recursive(f, a, c, TOL / 2, fa, fc, fd, approx_l) +
                    recursive(f, c, b, TOL / 2, fc, fb, fe, approx_r))

    fa = f(a); fb = f(b); fc = f((a + b) / 2)

    approx_0 = (b - a) * (fa + 4 * fc + fb) / 6

    return recursive(f, a, b, 10 * TOL, fa, fb, fc, approx_0)

def f_a(x):
    return 1 / (1 + x ** 4)

def f_a_i(t):
    return t ** 2 / (1 + t ** 4)

def f_b(x):
    return 1 / (1 + x ** 2) ** 3

def f_b_i(t):
    return t ** 4 / (t ** 2 + 1) ** 3

approx_a1 = adaptiveQuadrature(f_a, 0, 1, TOL)
approx_a2 = adaptiveQuadrature(f_a_i, 0, 1, TOL)
print(f"a: {approx_a1 + approx_a2}, {abs(approx_a1 + approx_a2 - math.pi / math.sqrt(2) / 2)}")

approx_b1 = adaptiveQuadrature(f_b, 0, 1, TOL)
approx_b2 = adaptiveQuadrature(f_b_i, 0, 1, TOL)
print(f"b: {approx_b1 + approx_b2}, {abs(approx_b1 + approx_b2 - 3 * math.pi / 16)}")