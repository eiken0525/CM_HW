end_point_1 = 1
end_point_2 = 4
tolerance = 0.001
maximum_number_of_iterations = 50
p_0 = 2.5
i = 1

def f(x):
    return x**3 + x - 4

FA = f(end_point_1)

while i <= maximum_number_of_iterations:
    p = end_point_1 + (end_point_2 - end_point_1) / 2
    FP = f(p)
    if (p == 0 or (end_point_2 - end_point_1) / 2 < tolerance):
        print(f"p = {p} with {i} iterations.")
        exit(0)
    i = i + 1
    if FA * FP > 0:
        end_point_1 = p
        FA = FP
    else:
        end_point_2 = p

print(f"Method failed after {maximum_number_of_iterations}.")