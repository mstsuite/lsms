import numpy as np
from io import StringIO

from scipy.integrate import simpson
from scipy.integrate import romb
from scipy.integrate import trapezoid

text = open('out.dat', "r").read()

s = StringIO(text)

D = np.genfromtxt(s)

r = D[:, 1]
y = D[:, 0]

res_1 = np.trapz(y, r)
res_2 = simpson(y, r)
res_4 = trapezoid(y, r)

print(res_1)
print(res_2)
print(res_4)
