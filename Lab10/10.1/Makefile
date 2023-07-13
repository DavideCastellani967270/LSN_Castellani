CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

Parallel_TSP.exe : Parallel_TSP.o random.o Parallel_TSP.h
	$(CC) random.o Parallel_TSP.o -o Parallel_TSP.exe
Parallel_TSP.o : Parallel_TSP.cpp Parallel_TSP.h
	$(CC) -c Parallel_TSP.cpp -o Parallel_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

run:
	mpiexec --hostfile hostfile -np 4 Parallel_TSP.exe

clean :
	rm *.o Parallel_TSP.exe *.out