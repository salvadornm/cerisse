import sympy as sp
import numpy as np
from fractions import Fraction


# Taylor series coefficents FD
AFD = np.array([
    [1, -3/2, 9/8, -27/48],
    [1, -1/2, 1/8, -1/48],
    [1,  1/2, 1/8,  1/48],
    [1,  3/2, 9/8,  27/48]
],  dtype=np.float64)

AFD_T = AFD.T

AFV = np.array([
    [1, -3/2, 7/6, -15/24],
    [1, -1/2, 1/6, -1/24],
    [1,  1/2, 1/6,  1/24],
    [1,  3/2, 7/6,  15/24]
], dtype=np.float64)

AFV_T = AFV.T

# Define vector B
B = np.array([0, 1, 0, 0], dtype=np.float64)

# Solve AT weight = B
weight = np.linalg.solve(AFD_T, B)

# Print solution
print("Solution weight (FINITE DIFFERENCES 4th ORDER):")
print(weight)

frac_arr = [Fraction(x).limit_denominator() for x in weight]

# Print the result
print("As fractions:")
print(frac_arr)

#-----------------------------------------------------------------
# Solve AT weight = B
weight_fv = np.linalg.solve(AFV_T, B)

# Print solution
print("Solution weight (FINITE VOLUME 4th ORDER):")
print(weight_fv)

frac_arr_fv = [Fraction(x).limit_denominator() for x in weight_fv]

# Print the result
print("As fractions:")
print(frac_arr_fv)

#------------------------------------------------------------------  6th order
# Taylor series coefficents FD
AFD = np.array([
    [1, -5/2, 25/8, -125/48, 625/384, -3125/3840],
    [1, -3/2,  9/8,  -27/48,  81/384,  -243/3840],
    [1, -1/2,  1/8,   -1/48,   1/384,    -1/3840],
    [1,  1/2,  1/8,    1/48,   1/384,     1/3840],
    [1,  3/2,  9/8,   27/48,  81/384,   243/3840],
    [1,  5/2, 25/8,  125/48, 625/384,  3125/3840]
],  dtype=np.float64)

AFD_T = AFD.T

# AFV = np.array([
#     [1, -3/2, 7/6, -15/24],
#     [1, -1/2, 1/6, -1/24],
#     [1,  1/2, 1/6,  1/24],
#     [1,  3/2, 7/6,  15/24]
# ], dtype=np.float64)

# AFV_T = AFV.T

# Define vector B
B = np.array([0, 1, 0, 0, 0, 0], dtype=np.float64)

# Solve AT weight = B
weight = np.linalg.solve(AFD_T, B)

# Print solution  
# [Fraction(-3, 640), Fraction(25, 384), Fraction(-75, 64), Fraction(75, 64), Fraction(-25, 384), Fraction(3, 640)]
print("Solution weight (FINITE DIFFERENCES 6th ORDER):")
print(weight)

frac_arr = [Fraction(x).limit_denominator() for x in weight]

# Print the result
print("As fractions:")
print(frac_arr)


