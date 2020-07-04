
esercizio10.2.exe : esercizio10.2.o random.o popolazione.o funzioni.o
	mpicxx random.o popolazione.o funzioni.o esercizio10.2.o -o esercizio10.2.exe 

esercizio10.2.o : esercizio10.2.cxx funzioni.h popolazione.h random.h
	mpicxx -c esercizio10.2.cxx -o esercizio10.2.o 

random.o : random.cpp random.h
	g++ -c random.cpp -o random.o

popolazione.o: popolazione.cxx popolazione.h random.h
	g++ -c -o popolazione.o popolazione.cxx

funzioni.o: funzioni.cxx funzioni.h random.h 
	g++ -c -o funzioni.o funzioni.cxx

