CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio3.1.exe : esercizio3.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio3.1.o -o esercizio3.1.exe
esercizio3.1.o : esercizio3.1.cxx funzioni.h
	$(CC) -c esercizio3.1.cxx -o esercizio3.1.o $(CFLAGS)
funzioni.o : funzioni.cxx funzioni.h
	$(CC) -c funzioni.cxx -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio3.1.exe seed.out


