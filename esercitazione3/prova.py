import matplotlib
import matplotlib.pyplot as plt
import numpy as np

i, costo, err = np.loadtxt("calldiretta.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
call=14.975790778311286
put=5.4595325819072364
plt.errorbar(i,costo,yerr=err)
plt.title("Costo opzione call")
plt.xlabel('numero step')
plt.ylabel('costo')
plt.grid(True)
N = 100
x = np.linspace(0.,100000.,N,endpoint=True)
y=np.zeros(N)
y=y+put
plt.plot(x,y,color="red") 
plt.show()
