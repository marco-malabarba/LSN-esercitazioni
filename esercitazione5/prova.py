import matplotlib
import matplotlib.pyplot as plt
import numpy as np

i, pi, err = np.loadtxt("media_100_gauss.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(i,pi,yerr=err)
plt.xlabel('numero tiri')
plt.ylabel('<r>')
plt.grid(True)
N = 100

#x = np.linspace(0.,10000000.,N,endpoint=True)
#y=np.zeros(N)
#y=y+np.pi
#plt.plot(x,y,color="red") 
plt.show()
