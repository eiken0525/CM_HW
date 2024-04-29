data = {1: 2.4142, 2: 2.6734, 3: 2.8974, 4: 3.0976, 5: 3.2804}

R_11 = 4 * (data[1]+data[5])
print(f"R_11 = {R_11}")

R_21 = 2 * (data[1]+data[3]+data[5])
print(f"R_21 = {R_21}")

R_31 = 1 * (data[1]+data[2]+data[3]+data[4]+data[5])
print(f"R_31 = {R_31}")

R_22 = R_21 + (R_21 - R_11) / 3
print(f"R_22 = {R_22}")

R_32 = R_31 + (R_31 - R_21) / 3
print(f"R_32 = {R_32}")

R_33 = R_32 + (R_32 - R_22) / 15
print(f"R_33 = {R_33}")