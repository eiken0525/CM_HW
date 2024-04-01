import math

def f(x):
    return - x**3 - math.cos(x)

p = [-1, 0]
q = [f(p[0]), f(p[1])]

for i in range(2):
    p_i = p[1]-f(p[1])*(p[1]-p[0])/(f(p[1])-f(p[0]))
    print(f"p({i+2}) = {p_i}")
    q_i = f(p_i)
    if q_i*q[1]<0:
        p[0] = p[1]
        q[0] = q[1]
    p[1] = p_i
    q[1] = q_i
