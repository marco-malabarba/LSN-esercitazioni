CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio1.1.exe : esercizio1.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio1.1.o -o esercizio1.1.exe
esercizio1.1.o : esercizio1.1.cxx
	$(CC) -c esercizio1.1.cxx -o esercizio1.1.o $(CFLAGS)
esercizio1.2.exe : esercizio1.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio1.2.o -o esercizio1.2.exe
esercizio1.2.o : esercizio1.2.cxx
	$(CC) -c esercizio1.2.cxx -o esercizio1.2.o $(CFLAGS)
esercizio1.3.exe : esercizio1.3.o funzioni.o random.o
	$(CC) random.o esercizio1.3.o funzioni.o -o esercizio1.3.exe
esercizio1.3.o : esercizio1.3.cxx
	$(CC) -c esercizio1.3.cxx -o esercizio1.3.o $(CFLAGS)
funzioni.o : funzioni.cxx funzioni.h
	$(CC) -c funzioni.cxx -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio1.1.exe seed.out


