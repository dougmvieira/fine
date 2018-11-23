from itertools import chain
from sympy import I, log, exp, sqrt, symbols
from sympy_utils import tree_constr, tree_codeblock, chain_rule_tree, diffs
from sympy.printing import fcode


fcode_opts = dict(standard=2008, source_format="free")

u = symbols('u')
t, kappa, a = symbols('t, kappa, a', positive=True)
rho = symbols('rho', real=True)
neg_t = symbols('-t', negative=True)
params = [t, kappa, a, rho]

beta = tree_constr('beta', (u, kappa, rho), kappa - rho*u*I)
d = tree_constr('d', (u, beta), sqrt(beta.symb**2 + u**2 + u*I))
g = tree_constr('g', (beta, d), (beta.symb - d.symb)/(beta.symb + d.symb))
h = tree_constr('h', (t, d, g), (g.symb*exp(-d.symb*t) - 1)/(g.symb - 1))
psi_1 = tree_constr('psi_1', (t, a, beta, d, h),
                    a*((beta.symb - d.symb)*t - 2*log(h.symb)))
psi_2 = tree_constr('psi_2', (t, beta, d, g),
                    (beta.symb - d.symb)*(1 - exp(-d.symb*t)
                                       )/(1 - g.symb*exp(-d.symb*t)))

psis = [psi_1, psi_2]
psi_grad = (chain_rule_tree(psi, diffs(psi), *params) for psi in psis)

with open("heston_psi.f90", "w") as f:
    print(fcode(tree_codeblock(*psis), **fcode_opts), file=f)
with open("heston_psi_grad.f90", "w") as f:
    print(fcode(tree_codeblock(*chain(psis, *psi_grad)), **fcode_opts), file=f)
