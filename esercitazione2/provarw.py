import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f(x,k):  
    return k*np.sqrt(x)

i, r, err = np.loadtxt("RWcontinuo.dat", usecols=(0,1,2), delimiter=' ', unpack='true')
plt.errorbar(i,r,yerr=err)
plt.title("Random walk discreto")
plt.xlabel('numero step')
plt.ylabel('$\sqrt{<|r_{step}|^2>_{RW}}$')
plt.grid(True)
p_opt, p_cov = curve_fit(f,i,r)
N = 100
x = np.linspace(1.,100.,N,endpoint=True)
y_fit = f(x,p_opt[0])

plt.plot(x,y_fit,color="red") 

print("Parametro k di best-fit:")
print(p_opt)

print("Incertezza sul parametro=")
print(np.sqrt(np.diagonal((p_cov))))

E_i=f(x,1)
sigma=err
chiquadro_i=np.zeros(np.size(x))
for i in range(np.size(x)):
	if sigma[i]==0:
		chiquadro_i[i]=0
	else: 
		chiquadro_i[i]=((E_i[i]-r[i])/sigma[i])**2
chiquadro=np.sum(chiquadro_i)
print(chiquadro)
plt.show()

