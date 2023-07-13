CC = g++
CFLAGS = -Wall -O3 --std=c++11

GA_Traveling_Salesman.exe : GA_Traveling_Salesman.o random.o
	$(CC) random.o GA_Traveling_Salesman.o -o GA_Traveling_Salesman.exe
GA_Traveling_Salesman.o : GA_Traveling_Salesman.cpp GA_Traveling_Salesman.h random.h
	$(CC) -c GA_Traveling_Salesman.cpp -o GA_Traveling_Salesman.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GA_Traveling_Salesman.exe *.out
run:
	./GA_Traveling_Salesman.exe