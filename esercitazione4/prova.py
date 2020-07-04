import matplotlib
import matplotlib.pyplot as plt
import numpy as np

misura = np.loadtxt("output_temp.dat", usecols=(0), delimiter=' ', unpack='true')
j=np.arange(len(misura))
j+=1

plt.plot(j,misura)

plt.title('Andamento temperatura')
plt.xlabel('Misura')
plt.ylabel('$T^*$')
plt.grid(True)
plt.show()
