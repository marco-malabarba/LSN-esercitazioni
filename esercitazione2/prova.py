import matplotlib
import matplotlib.pyplot as plt
import numpy as np

i, pi, err = np.loadtxt("integrale3.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(i,pi,yerr=err)
plt.xlabel('numero tiri')
plt.ylabel('I')
plt.grid(True)
N = 100
x = np.linspace(0.,1000000.,N,endpoint=True)
y=np.zeros(N)
y=y+1.
plt.plot(x,y,color="red") 
plt.show()
