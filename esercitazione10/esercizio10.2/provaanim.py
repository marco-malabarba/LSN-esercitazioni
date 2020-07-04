import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def update(i):
    xdata.append(x[i])
    ydata.append(y[i])
    ln.set_data(xdata, ydata)
    return ln,

fig, ax = plt.subplots()
x,y = np.loadtxt("percorso_quad.dat", usecols=(0,1), delimiter=' ', unpack='true')
ax.set_xlim(-12.5,12.5)
ax.set_ylim(-12.5,12.5)
plt.title("Percorso")
plt.plot(x,y,'ro')

xdata, ydata = [], []
ln, = plt.plot([], [], marker='.',markerfacecolor="r")
ani = FuncAnimation(fig, update, frames= len(x), blit=False, interval=1000, repeat=False)

plt.grid(True)
plt.show()
