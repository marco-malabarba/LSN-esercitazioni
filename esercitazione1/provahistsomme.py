import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(x,var,mu,k):  
    return k/(np.sqrt(2*np.pi*var)) * np.exp(-((x-mu)**2)/(2*var))

data = np.loadtxt("somme.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 20, range=(2.5,3.5),color='red',ec='grey')
data_entries, bins = np.histogram(data, bins=bins)
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
N = 100
x = np.linspace(2.5,3.5,N,endpoint=True)
p_opt, p_cov = curve_fit(f, xdata=binscenters, ydata=data_entries)
y_fit = f(x,p_opt[0],p_opt[1],p_opt[2])
plt.plot(x,y_fit) 
plt.title('Distribuzione $N=200$')
plt.xlabel('$S_{200}$')
plt.ylabel('occorrenze')
plt.grid(True)
plt.show()


print("Parametri di best-fit per la Gaussiana [varianza,mu,k]=")
print(p_opt)
print("Incertezza sui parametri=")
print(np.sqrt(np.diagonal(p_cov)))

E_i=f(binscenters,p_opt[0],p_opt[1],p_opt[2])
sigma_chiquadro=np.sqrt(E_i)
chiquadro_i=((data_entries-E_i)/(sigma_chiquadro))**2

chiquadro=np.sum(chiquadro_i)
print("chi-quadro:")
print(chiquadro)

