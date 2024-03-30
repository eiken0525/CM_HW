import math

def adm(p1, p2, p3):
    return p1 - (p2 - p1) ** 2 / (p3 - 2 * p2 + p1)

def f(x):
    return math.cos(x)

p = [0.5]
p_hat = []

for i in range(7):
    p.append(f(p[-1]))

for i in range(5):
    p_hat.append(adm(p[i], p[i+1], p[i+2]))

print(p_hat)