import math

a = 0
b = math.pi
n = 1000

h = (b - a) / n

def x(n):
    return a + n * h

def f(x):
    a = math.sqrt(4 + 5 * math.sin(x) ** 2)
    return a

sum = f(a) + f(b)

for i in range(n//2 - 1):
    sum += 2 * f(x(2 * i))

for i in range(n//2):
    sum += 4 * f(x(2 * i-1))

print(2 * h * sum / 3)



