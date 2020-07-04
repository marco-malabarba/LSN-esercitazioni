import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def f(x,mu,sigma): 
	psi=np.exp(-(x-mu)**2/(2*(sigma**2))) + np.exp(-(x+mu)**2/(2*(sigma**2))) #funzione d'onda non normalizzata
	norm=2*np.sqrt(np.pi)*sigma*np.exp(-mu*mu/(sigma*sigma))+2*np.sqrt(np.pi)*sigma; #fattore di normalizzazione
	return (psi**2)/norm


i, pi, err = np.loadtxt("ground_state.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(i,pi,yerr=err)
plt.xlabel('numero blocchi')
plt.ylabel('$E_{gs}$')
plt.grid(True)
plt.show()

mu=0.802393
sigma=0.618400
N = 1000
x = np.linspace(-3,3,N,endpoint=True)
fx=f(x,mu,sigma)
x_med,occ=np.loadtxt("histo.dat", usecols=(0,1), delimiter=' ', unpack='true')
plt.hist(x_med, weights=occ,bins=len(x_med), density=1)
plt.plot(x,fx)

def Vpot(x):
    return (x**2 - 2.5)*x**2
    #return 0.5*x**2

hbar = 1
m = 1
a = 10
N = 1000 # number of iterations

# Step sizes
x = np.linspace(-a/2, a/2, N)
dx = x[1] - x[0] # the step size
V = Vpot(x)

# The central differences method: f" = (f_1 - 2*f_0 + f_-1)/dx^2

CDiff = np.diag(np.ones(N-1),-1)-2*np.diag(np.ones(N),0)+np.diag(np.ones(N-1),1)
# np.diag(np.array,k) construct a "diagonal" matrix using the np.array
# The default is k=0. Use k>0 for diagonals above the main diagonal, 
# and k<0 for diagonals below the main diagonal

# Hamiltonian matrix
H = (-(hbar**2)*CDiff)/(2*m*dx**2) + np.diag(V)

# Compute eigenvectors and their eigenvalues
E,psi = np.linalg.eigh(H)

# Take the transpose & normalize
psi = np.transpose(psi)
psi = psi/np.sqrt(dx)
plt.plot(x,(psi[0])**2)



plt.grid(True)
plt.show()
