CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio9.1.exe : esercizio9.1.o random.o popolazione.o funzioni.o
	$(CC) random.o popolazione.o funzioni.o esercizio9.1.o -o esercizio9.1.exe

esercizio9.1.o : esercizio9.1.cxx funzioni.h popolazione.h random.h
	$(CC) -c esercizio9.1.cxx -o esercizio9.1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

popolazione.o: popolazione.cxx popolazione.h random.h
	g++ -c -o popolazione.o popolazione.cxx

funzioni.o: funzioni.cxx funzioni.h random.h 
	g++ -c -o funzioni.o funzioni.cxx

clean :
	rm *.o esercizio9.1.exe seed.out
