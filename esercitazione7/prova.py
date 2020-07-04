import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math



r,gr,err = np.loadtxt("output.gave.gas.0", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(r,gr,yerr=err, color ="red")

r,gr,err = np.loadtxt("output.gave_md.0", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(r,gr,yerr=err, color ="blue")

plt.grid(True)
plt.show()
