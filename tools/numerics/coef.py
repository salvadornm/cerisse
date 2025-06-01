import sympy as sp

# Define symbols
h = sp.Symbol('h')
x = sp.Symbol('x')
f = sp.Function('f')

# Expansion point (face-centered): x_{i+1/2}
x0 = sp.Symbol('x_{i+1/2}')  # for LaTeX clarity

# Relative offsets from x_{i+1/2}
offsets = {
    'i+3': 2.5*h,
    'i+2': 1.5*h,
    'i+1': 0.5*h,
    'i'  : -0.5*h,
    'i-1': -1.5*h,
    'i-2': -2.5*h
}

# Order of expansion
order = 5

# Expand and collect LaTeX output
for label, delta in offsets.items():
    taylor = f(x0)
    for n in range(1, order + 1):
        term = (delta**n / sp.factorial(n)) * sp.diff(f(x), x, n)
        taylor += term
    taylor = sp.simplify(taylor)
    latex_expr = sp.latex(sp.Eq(f(sp.Symbol(label)), taylor))
    print(f"\\[{latex_expr}\\]\n")

