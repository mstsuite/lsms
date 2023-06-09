import matplotlib.pyplot as plt
import numpy as np

energies = []

for i in range(7):
    res = np.loadtxt(f"{i}/k.out")
    energies.append(res[-1, 1])

lat = np.array(
    [
        3.75,
        3.85,
        3.90,
        3.95,
        4.00,
        4.05,
        4.15
    ]
)

plt.plot(lat, energies, ls='--', marker='x', color="red")
plt.savefig("res.png")
