import math

def f(x):
    return - x**3 - math.cos(x)

p = [-1, 0]
for i in range(2):
    p.append(p[i+1]-f(p[i+1])*(p[i+1]-p[i])/(f(p[i+1])-f(p[i])))
    print(f"p({i+2}) = {p[-1]}")