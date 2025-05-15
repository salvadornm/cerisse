import sympy as sp

# Symbols
x, t = sp.symbols('x t')
gamma = sp.Rational(7, 5)  # gamma = 1.4

# Manufactured solution
# time dependent
#rho = 1 + 0.2 * sp.sin(2*sp.pi * x) * sp.cos(sp.pi * t)
#u = 1 + 0.1 * sp.cos(2*sp.pi * x) * sp.sin(sp.pi * t)
#p = 1 + 0.3 * sp.cos(2*sp.pi * x) * sp.cos(sp.pi * t)
# steady
rho = 1 + 0.2 * sp.sin(2*sp.pi * x) 
u = 1
p = 1 + 0.3 * sp.cos(2*sp.pi * x)

# Derived quantities
rho_u = rho * u
E = p / (gamma - 1) + 0.5 * rho * u**2

# Fluxes
F1 = rho_u
F2 = rho * u**2 + p
F3 = u * (E + p)

# Time derivatives
dU1_dt = sp.diff(rho, t)
dU2_dt = sp.diff(rho_u, t)
dU3_dt = sp.diff(E, t)

# Spatial derivatives of fluxes
dF1_dx = sp.diff(F1, x)
dF2_dx = sp.diff(F2, x)
dF3_dx = sp.diff(F3, x)

# Forcing terms: S = ∂U/∂t + ∂F/∂x
S1 = dU1_dt + dF1_dx
S2 = dU2_dt + dF2_dx
S3 = dU3_dt + dF3_dx

# Simplify
S1_simplified = sp.simplify(S1)
S2_simplified = sp.simplify(S2)
S3_simplified = sp.simplify(S3)

# Output
print("rho :")
print(sp.latex(rho))

print("u :")
print(sp.latex(u))

print("p :")
print(sp.latex(p))


# Output forcing terms
print("S_rho (mass equation: \n")
print(sp.latex(S1))

print("----------------------- \n")

print("S_momentum (momentum equation): \n")
#sp.pprint(S2_simplified)
print(sp.latex(S2))

print("-----------------------\n")

print("S_energy (energy equation):  \n")
#sp.pprint(S3_simplified)
print(sp.latex(S3))


