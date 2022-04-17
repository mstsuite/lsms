import numpy as np
from scipy.integrate import simps

a = np.loadtxt("libxc/example.txt")

ir_max = 501

val = a[:ir_max,0]
r = a[:ir_max,1]

a1 = np.trapz(val, r)

print(a1)

a2 = simps(val, r)

print(a2)
