import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math


def grafici(stato):
	intestazione="Stato " + stato
	print(intestazione)
	
	if stato=="solido":
		titolo1="sol"
	if stato=="liquido":
		titolo1="liq"
	if stato=="gassoso":
		titolo1="gas"
		
	nome_file_1="output.epot."+titolo1+".0"
	
	bl,energia,err = np.loadtxt(nome_file_1, usecols=(0,2,3), delimiter=' ', unpack='true')
	
	energia=energia*120.*1.380649*(10**(-23)) #conversione in SI
	err=err*120.*1.380649*(10**(-23))
	
	plt.errorbar(bl,energia,yerr=err)
	plt.title("Andamento energia potenziale per particella")
	plt.xlabel("numero blocchi")
	plt.ylabel("U/N [J]")
	plt.grid(True)
	plt.show()
	
	nome_file_2="output.pres."+titolo1+".0"
	bl,pres,err = np.loadtxt(nome_file_2, usecols=(0,2,3), delimiter=' ', unpack='true')
	
	pres=pres*120.*1.380649*(10**(-23))/((0.34*(10**(-9)))**3)
	err=err*120.*1.380649*(10**(-23))/((0.34*(10**(-9)))**3)
	
	plt.errorbar(bl,pres,yerr=err)	
	plt.title("Andamento pressione")
	plt.xlabel("numero blocchi")
	plt.ylabel("P[Pa]")
	plt.grid(True)
	plt.show()	
	
grafici("solido")
grafici("liquido")
grafici("gassoso")
	
	
