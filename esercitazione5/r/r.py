from sympy import *
from sympy.abc import x,y,z,r,theta
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math


def fattoriale(n):
	k=1
	while n>1:
		k=k*n
		n=n-1
	return k


n=int(input("inserire il numero quantico principale: "))
l=int(input("inserire il numero quantico azimutale: "))
ml=int(input("inserire il numero quantico magnetico: "))

filename="nlm.dat"
myfile = open(filename, 'w')
myfile.write(str (n))
myfile.write(' ')
myfile.write(str (l))
myfile.write(' ')
myfile.write(str (ml))
myfile.close()

#trascuro tutte le normalizzazioni (in Metropolis si elidono)

phi=1 #trascuro il fattore exp(i...): nel modulo quadro sparisce

i=0
Pl=((x**2)-1.)**l

while i<l:
	Pl=diff(Pl,x)
	i+=1

i=0

while i<abs(ml):
	Pl=diff(Pl,x)
	i+=1

Pl=Pl*(1-(x)**2)**(abs(ml)/2.)
Plm=Pl.subs({x:cos(theta)})


THETA=Plm

i=0
Lp=(x**(n+l))*exp(-x)
while i<n+l:
	Lp=diff(Lp,x)
	i=i+1

Lp=Lp*exp(x)

simplify(Lp)
expand(Lp)
powsimp(Lp)


i=0
while i<2*l+1:
	Lp=diff(Lp,x)
	i=i+1

simplify(Lp)
expand(Lp)
powsimp(Lp)

i=0
Lpq=Lp.subs({x:(2.*r/n)})

simplify(Lpq)
expand(Lpq)
powsimp(Lpq)

R=(((2/n)*r)**l)*exp(-r/n)*Lpq
simplify(R)
expand(R)
powsimp(R)



PSI=(R*THETA*phi)**2 #modulo quadro funzione d'onda in coordinate sferiche
psi=PSI.subs({r:sqrt(x*x+y*y+z*z), theta: acos(z/sqrt(x*x+y*y+z*z))}) #modulo quadro funzione d'onda in coordinate cartesiane

simplify(psi)
expand(psi)

filename = "funzione.cxx"
myfile = open(filename, 'w')
myfile.write("#include <iostream> \n")
myfile.write("#include <cmath> \n")
myfile.write("using namespace std;")

myfile.write("double psi(double x, double y, double z) {\n")
myfile.write("return " + ccode(psi, standard='C89') + ";\n" )
myfile.write("}")
myfile.close()







