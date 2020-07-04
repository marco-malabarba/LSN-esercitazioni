import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x,y = np.loadtxt("percorso_0.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_1.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_2.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()

x,y = np.loadtxt("percorso_3.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.plot(x,y,marker='.',markerfacecolor="r")
plt.title("Percorso")
plt.grid(True)
plt.show()
