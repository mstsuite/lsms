import numpy as np

a = np.loadtxt("built_in/outfile6")
b = np.loadtxt("libxc/outfile6")

a1 = np.trapz(a)
b1 = np.trapz(b)

print(a1)
print(b1)
