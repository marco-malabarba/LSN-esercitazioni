CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio2.1.exe : esercizio2.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio2.1.o -o esercizio2.1.exe
esercizio2.1.o : esercizio2.1.cxx funzioni.h
	$(CC) -c esercizio2.1.cxx -o esercizio2.1.o $(CFLAGS)
esercizio2.2.exe : esercizio2.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio2.2.o -o esercizio2.2.exe
esercizio2.2.o : esercizio2.2.cxx funzioni.h
	$(CC) -c esercizio2.2.cxx -o esercizio2.2.o $(CFLAGS)
funzioni.o : funzioni.cxx funzioni.h
	$(CC) -c funzioni.cxx -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio2.1.exe seed.out


