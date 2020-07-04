LSN\_exercises\_delivery

Sono presenti 12 cartelle, una per ogni esercitazione. All'interno di ciascuna vi sono: il notebook con il testo degli esercizi; il notebook in cui vengono commentati i vari risultati (intitolati "esercitazioneXY.ipynb", con XY=numero esercitazione); i programmi con cui si sono ottenuti i risultati; i file di output. 
Tutti i programmi in python, generalmente intitolati "prova.py", sono analoghi (o molto simili) a quelli presenti negli specchietti dei jupyter-notebook.

# ISTRUZIONI PER COMPILARE

## Esercitazione 1

I tre esercizi sono risolti con i programmi in C++ "esercizio1.1.cxx", "esercizio1.2.cxx" e "esercizio1.3.cxx". Nella cartella è presente un makefile. Per compilare "esercizio1.1.cxx" bisogna utilizzare il comando "make esercizio1.1.exe". Per gli altri due il comando è analogo, è sufficiente modificare il numero dell'esercizio.

## Esercitazioni 2, 3, 5 
 
Valgono le stesse considerazioni dell'esercitazione 1
 
## Esercitazione 4

Il programma è intitolato "MolDyn\_NVE.cpp". Per compilare utilizzare la classica istruzione "g++ MolDyn\_NVE.cpp -o MolDyn\_NVE.exe". Per iniziare la simulazione con i parametri desiderati modificare il file "input.dat".

## Esercitazione 6

Il programma è intitolato "Monte\_Carlo\_ISING\_1D.cpp". Per compilare è sufficiente scrivere a terminale "make": viene creato l'eseguibile "ising.exe". Per iniziare la simulazione con i parametri desiderati, modificare il file "input.dat".

## Esercitazione 7

Il programma è intitolato "Monte\_Carlo\_NVT.cpp". Per compilare è sufficiente scrivere a terminale "make": viene creato l'eseguibile "Monte\_Carlo\_NVT.exe". Per iniziare la simulazione con i parametri desiderati, modificare il file "input.dat". Nella cartella è presente anche il programma "MolDyn\_NVE.cpp" dell'esercitazione 4, necessario per risolvere gli esercizi 7.3 e 7.4.

## Esercitazione 8

I primi due esercizi sono risolti con il programma "variazionale.cxx". Per creare l'eseguibile "variazionale.exe" utilizzare il comando "g++ variazionale.cxx -o variazionale.exe". Per il terzo programma è sufficiente, invece, l'istruzione "make" e, per iniziare la simulazione con i parametri desiderati, modificare il file "input.dat". 

## Esercitazione 9

La classe è intitolata "popolazione", il main del programma "esercizio9.1.cxx". Con l'istruzione "make" viene creato l'eseguibile "esercizio9.1.exe"

## Esercitazione 10

Per ciascuno dei due esercizi si può trovare una sottocartella. All'interno di entrambe è presente un makefile. Per compilare è sufficiente, quindi, "make". I nomi degli eseguibili sono "esercizio10.1.exe" e "esercizio10.2.exe". Per l'esercizio 10.2 si sono utilizzate le librerie MPI, quindi per lanciare il programma è necessaria l'istruzione "mpiexec -np 4 esercizio10.2.exe". 

 
 
