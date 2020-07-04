import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x,y = np.loadtxt("percorso_0_nomig.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_1_nomig.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_2_nomig.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_3_nomig.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()
