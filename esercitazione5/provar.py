import matplotlib
import matplotlib.pyplot as plt
import numpy as np

r=np.loadtxt("andamento.dat", usecols=(0), delimiter=' ', unpack='true')
rfar=np.loadtxt("andamento_far.dat", usecols=(0), delimiter=' ', unpack='true')
i=np.arange(len(r))
i=i+1
plt.title('andamento r')
plt.xlabel('numero step')
plt.ylabel('r')

plt.plot(i,r,color="red",label="partenza da un punto \'ragionevole\'")
plt.plot(i,rfar,color="blue",label="partenza molto lontano dall'origine")
plt.grid(True)
legend = plt.legend(loc='upper center',shadow=True)
plt.show()

