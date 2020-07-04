import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def grafici(stato):
	if stato=="solido":
		titolo1="solido"
	if stato=="liquido":
		titolo1="liq"
	if stato=="gassoso":
		titolo1="gas"
		
	titolo2="temp"
	titolo="average_"+titolo2+"_"+titolo1+".dat"
	
	step, t,err = np.loadtxt(titolo, usecols=(0,1,2), delimiter=' ', unpack='true')
	t=t*120.
	err=err*120.
	plt.errorbar(step,t,yerr=err)
	plt.title ('Temperatura')
	plt.xlabel('numero di misure')
	plt.ylabel('T[K]')
	plt.grid(True)
	plt.show()
	
	titolo2="ekin"
	titolo="average_"+titolo2+"_"+titolo1+".dat"
		
	step, ekin,err = np.loadtxt(titolo, usecols=(0,1,2), delimiter=' ', unpack='true')
	ekin=ekin*120.*1.380649*(10**(-23))
	err=err*120.*1.380649*(10**(-23))
	plt.errorbar(step,ekin,yerr=err)
	plt.title ('Energia cinetica per particella')
	plt.xlabel('numero di misure')
	plt.ylabel('$K/N[J]$')
	plt.grid(True)
	plt.show()
	
	titolo2="epot"
	titolo="average_"+titolo2+"_"+titolo1+".dat"
	
	step, epot,err = np.loadtxt(titolo, usecols=(0,1,2), delimiter=' ', unpack='true')
	epot=epot*120.*1.380649*(10**(-23))
	err=err*120.*1.380649*(10**(-23))
	plt.errorbar(step,epot,yerr=err)
	plt.title ('Energia potenziale per particella')
	plt.xlabel('numero di misure')
	plt.ylabel('U/N [J]')
	plt.grid(True)
	plt.show()
	
	titolo2="etot"
	titolo="average_"+titolo2+"_"+titolo1+".dat"
	
	step, etot,err = np.loadtxt(titolo, usecols=(0,1,2), delimiter=' ', unpack='true')
	etot=etot*120.*1.380649*(10**(-23))
	err=err*120.*1.380649*(10**(-23))
	plt.errorbar(step,etot,yerr=err)
	plt.title ('Energia totale per particella')
	plt.xlabel('numero di misure')
	plt.ylabel('E/N [J]')
	plt.ticklabel_format(useOffset=False)
	plt.grid(True)
	plt.show()
	

grafici("gassoso")
