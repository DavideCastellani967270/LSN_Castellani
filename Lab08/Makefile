CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_VARIATIONAL.exe : Monte_Carlo_VARIATIONAL.o random.o
	$(CC) random.o Monte_Carlo_VARIATIONAL.o -o Monte_Carlo_VARIATIONAL.exe
Monte_Carlo_VARIATIONAL.o : Monte_Carlo_VARIATIONAL.cpp
	$(CC) -c Monte_Carlo_VARIATIONAL.cpp -o Monte_Carlo_VARIATIONAL.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_VARIATIONAL.exe *.out
run:
	./Monte_Carlo_VARIATIONAL.exe