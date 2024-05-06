N1 = [1.570796, 1.896119, 1.974232, 1.993570]
N2 = []
N3 = []
N4 = []

for i in range(3):
    N2.append(N1[i+1] + (N1[i+1] - N1[i]) / 3)

for i in range(2):
    N3.append(N2[i+1] + (N2[i+1] - N2[i]) / 15)

for i in range(1):
    N4.append(N3[i+1] + (N3[i+1] - N3[i]) / 63)

print(N1)
print(N2)
print(N3)
print(N4)