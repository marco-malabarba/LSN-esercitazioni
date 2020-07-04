CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio5.1.exe : esercizio5.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio5.1.o -o esercizio5.1.exe
esercizio5.1.o : esercizio5.1.cxx funzioni.h
	$(CC) -c esercizio5.1.cxx -o esercizio5.1.o $(CFLAGS)
esercizio1.2.exe : esercizio1.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio1.2.o -o esercizio1.2.exe
funzioni.o : funzioni.cxx funzioni.h
	$(CC) -c funzioni.cxx -o funzioni.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio5.1.exe seed.out


