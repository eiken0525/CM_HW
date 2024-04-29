import math

a = 0
b = 2
n = 20000

h = (b-a)/n

def x(i):
    return a + h*i

def f(x):
    return math.e ** (2 * x) * math.sin(3 * x)

ctSum = f(a)+f(b)
for i in range(n-1):
    ctSum += 2 * f(x(i+1))
print("Composite trapezoidal rule:", h * ctSum / 2)

cSSum = f(a)+f(b)
for i in range(n//2-1):
    cSSum += 2 * f(x(2*(i+1)))
for i in range(n//2):
    cSSum += 4 * f(x(2*(i+1)-1))
print("Composite Simpson's rule:", h * cSSum / 3)


def x_cm(i):
    return a + h*(i+1)
n = 19998
h = (b-a)/(n+2)
cmSum = 0
for i in range(n//2+1):
    cmSum += f(x_cm(2*i))
print("Composite midpoint rule:", 2 * h * cmSum)


