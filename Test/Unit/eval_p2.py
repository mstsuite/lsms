import numpy as np
from scipy.integrate import simps

for i in range(1,4):

    a = np.loadtxt(f"ezrho{i}.out")

    ir_max = 501
    val = a[:ir_max,0]
    r = a[:ir_max,1]
    a2 = simps(val, r)
    print(a2)

