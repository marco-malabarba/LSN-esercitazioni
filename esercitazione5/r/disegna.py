import matplotlib
import matplotlib.pyplot as plt
import numpy as np

i, r, err = np.loadtxt("risultati.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(i,r,yerr=err)
plt.xlabel('numero step')
plt.ylabel('<r>/$a_0$')
plt.grid(True)

n,l,ml=np.loadtxt("nlm.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
N = 100
x = np.linspace(0.,10000000.,N,endpoint=True)
y=np.zeros(N)
y=y+(3*n*n-l*(l+1))/2.
plt.plot(x,y,color="red") 

plt.show()
