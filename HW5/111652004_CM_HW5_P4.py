def eulerMethod(f, h, a, b):
    y_0 = 1
    N = int((b - a) / h)
    for i in range(N):
        y_0 += h * f(a + h * (i + 1), y_0)
    return y_0

def f_4(t, y):
    return -y + t + 1

for h in [0.2, 0.1, 0.05]:
    print(f"With h = {h}, the approximation is {eulerMethod(f_4, h, 0, 5)}.")


    