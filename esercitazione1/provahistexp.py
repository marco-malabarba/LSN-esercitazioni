import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(x,var,mu,k):  
    return k/(np.sqrt(2*np.pi*var)) * np.exp(-((x-mu)**2)/(2*var))



data = np.loadtxt("exp1.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 12, range=(0,6),ec='grey')
plt.title('Distribuzione esponenziale $N=1$')
plt.xlabel('$S_1$')
plt.ylabel('occorrenze')
plt.grid(True)
plt.show()
data = np.loadtxt("exp2.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 12, range=(0,6),color='green',ec='grey')
plt.title('Distribuzione esponenziale $N=2$')
plt.xlabel('$S_2$')
plt.ylabel('occorrenze')
plt.grid(True)
plt.show()
data = np.loadtxt("exp10.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 15, range=(0,2.5),color='purple',ec='grey')
plt.title('Distribuzione esponenziale $N=10$')
plt.xlabel('$S_{10}$')
plt.ylabel('occorrenze')
plt.grid(True)
plt.show()

data = np.loadtxt("exp100.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 20, range=(0.5,1.5),color='red',ec='grey')
data_entries, bins = np.histogram(data, bins=bins)
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
N = 100
x = np.linspace(0.5,1.5,N,endpoint=True)
p_opt, p_cov = curve_fit(f, xdata=binscenters, ydata=data_entries)
y_fit = f(x,p_opt[0],p_opt[1],p_opt[2])
plt.plot(x,y_fit) 

plt.title('Distribuzione esponenziale $N=100$')
plt.xlabel('$S_{100}$')
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
print("chiquadro:")
print(chiquadro)




data = np.loadtxt("exp100.dat", usecols=(0), delimiter=' ', unpack='true')
n, bins, patches = plt.hist(data, 20, range=(0.8,1.2),color='red',ec='grey')
data_entries, bins = np.histogram(data, bins=bins)
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
N = 100
x = np.linspace(0.8,1.2,N,endpoint=True)
p_opt, p_cov = curve_fit(f, xdata=binscenters, ydata=data_entries)
y_fit = f(x,p_opt[0],p_opt[1],p_opt[2])
plt.plot(x,y_fit)
plt.title('Distribuzione esponenziale $N=100$') 
plt.xlabel('$S_{100}$')
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
print("chiquadro:")
print(chiquadro)


