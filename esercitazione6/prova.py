import matplotlib
import matplotlib.pyplot as plt
import numpy as np

misura = np.loadtxt("output.en_ist.0", usecols=(0), delimiter=' ', unpack='true')
j=np.arange(len(misura))
j+=1

plt.plot(j,misura)

plt.title('Andamento energia')
plt.xlabel('Misura')
plt.ylabel('$H/N$')
plt.grid(True)
plt.show()
