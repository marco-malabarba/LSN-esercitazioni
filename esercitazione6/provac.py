import matplotlib
import matplotlib.pyplot as plt
import numpy as np

t, mis, err= np.loadtxt("ene_metro.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="red",marker='*',label="Metropolis")
t, mis, err= np.loadtxt("ene_gibbs.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="blue",marker='*',label="Gibbs")
plt.title('andamento U/N(T)')
plt.xlabel('U/N')
plt.ylabel('T')
points=500
T = np.linspace(0.2,3.0,num=points)
beta = 1/T
J = 1.0
Ns = 50
th = np.tanh(J/T)
thN= th**Ns
ch = 1/th
e = -J*( th + ch*thN )/( 1 + thN )
plt.plot(T, e,color="black",label="curva teorica")

plt.grid(True)
legend = plt.legend(loc='upper center',shadow=True)

plt.show()

t, mis, err= np.loadtxt("heat_metro.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="red",marker='*',label="Metropolis")
t, mis, err= np.loadtxt("heat_gibbs.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="blue", marker='*',label="Gibbs")
plt.title('andamento C(T)')
plt.xlabel('C')
plt.ylabel('T')
heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)
plt.plot(T, heat,color="black",label="curva teorica")

plt.grid(True)
legend = plt.legend(loc='upper center',shadow=True)

plt.show()

t, mis, err= np.loadtxt("mag_metro.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="red",marker='*',label="Metropolis")
t, mis, err= np.loadtxt("mag_gibbs.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="blue",marker='*',label="Gibbs")
plt.title('andamento M(T)')
plt.xlabel('T')
plt.ylabel('M')

h=0.02
b = 1/T

l1 = np.exp(b*J)*np.cosh(b*h)+np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
l2 = np.exp(b*J)*np.cosh(b*h)-np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
Z = l1**Ns + l2**Ns
M = (np.exp(b*J)*np.sinh(b*h)*((l1**(Ns-1))*(1+np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))) 
        + (l2**(Ns-1))*(1-np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J)))))/(Z)
plt.plot(T, M,color="black",label="curva teorica")
plt.grid(True)
legend = plt.legend(loc='upper center',shadow=True)

plt.show()

t, mis, err= np.loadtxt("chi_metro.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="red",marker='*',label="Metropolis")
t, mis, err= np.loadtxt("chi_gibbs.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(t,mis,yerr=err,color="blue",marker='*',label="Gibbs")
plt.title('andamento $\chi$(T)')
plt.xlabel('T')
plt.ylabel('$\chi$')

X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
plt.plot(T, X, color="black", label="curva teorica")
plt.grid(True)
legend = plt.legend(loc='upper center',shadow=True)

plt.show()


