p_0 = 0.6
tolerance = 0.0001
maximum_number_of_iterations = 50
i = 1

def g(x):
    return 0.5**x

while i <= maximum_number_of_iterations:
    p = g(p_0)
    if abs(p - p_0) < tolerance:
        print(f"p = {p} with {i} iterations.")
        exit(0)
    i = i + 1
    p_0 = p

print(f"Method failed after {maximum_number_of_iterations}.")