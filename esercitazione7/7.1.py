import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math


def error(AV,AV2,n):  # Function for statistical uncertainty estimation
	if n==0:
		return 0
	else:
		return math.sqrt((AV2 - AV**2)/n)

def function(stato,grandezza):
	intestazione1="Stato "+stato
	intestazione2=grandezza
	print(intestazione1)
	print(intestazione2)
	if stato=="solido":
		titolo1="sol"
	if stato=="liquido":
		titolo1="liq"
	if stato=="gassoso":
		titolo1="gas"
	
	titolo="istantanei_"+titolo1+"_"+grandezza+".dat"
	
	energia = np.loadtxt(titolo, usecols=(0), delimiter=' ', unpack='true') 
	energia_2=energia**2

	media=np.sum(energia)/len(energia)
	sigma=np.sum(energia_2)/len(energia)
	sigma=sigma-media**2



	autocorr=np.zeros(200)


	for i in range(200):
		counter=0
		for j in range(len(energia)-i):
			autocorr[i]+=energia[j]*energia[i+j]
			counter=counter+1
		autocorr[i]=autocorr[i]/(counter)
		autocorr[i]=(autocorr[i]-media**2)/sigma
		#print(autocorr[i])


	x=np.arange(200)	
	plt.plot(x,autocorr)
	plt.grid(True)
	plt.title("Andamento autocorrelazione")
	plt.xlabel("$\tau$")
	plt.ylabel("autocorrelazione")
	plt.show()

	print("valori con errori")
	L=[10,20,50,100,200,500,1000,2000,5000]

	for z in range(len(L)):
		
		N=int(len(energia)/L[z])
		ave=np.zeros(N)
		av2=np.zeros(N)
		for i in range(N):
			sum = 0
			for j in range(L[z]):
				k = j+i*L[z]
				sum += energia[k]
			ave[i] = sum/L[z]      
			av2[i] = (ave[i])**2 
	
		sum_prog=0
		su2_prog=0
		for i in range(N):
			sum_prog += ave[i] 
			su2_prog += av2[i]
		sum_prog/=(N) 
		su2_prog/=(N) 
		err_prog = error (sum_prog,su2_prog,N+1)
		#print(sum_prog," ",err_prog)
		plt.errorbar(L[z],sum_prog,yerr=err_prog,color="black")
	
	plt.xscale('log')
	title="Andamento errori sulla "+grandezza
	plt.title(title)
	plt.xlabel("L")
	plt.ylabel(grandezza)
	plt.grid(True)
	plt.show()


function("solido","energia")
function("solido","pressione")
function("liquido","energia")
function("liquido","pressione")
function("gassoso","energia")
function("gassoso","pressione")

