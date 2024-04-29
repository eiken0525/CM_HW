import math

def f(t):
    return (6 * t + 26) / (3 * t ** 2 + 26 * t - 77)

n_values = []

node_data = {
    2: ([0.57735, -0.57735], [1, 1]),
    3: ([0.77460, 0, -0.77460], [0.55556, 0.88889, 0.55556]),
    4: ([0.86114, 0.33998, -0.33998, -0.86114], [0.34785, 0.65215, 0.65215, 0.34785])
}

for num_nodes in range(2, 5):
    sum = 0
    roots, coeffs = node_data[num_nodes]
    for i in range(num_nodes):
        sum += coeffs[i] * f(roots[i])
    n_values.append(sum)

rel_error = []

for i in range(3):
    rel_error.append(abs(n_values[i] - math.log(0.48)) / abs(math.log(0.48)))

print(n_values)
print(rel_error)

