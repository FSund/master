#CC = icpc -I../mts0_io -O3
CC = g++ -I../mts0_io -std=c++11

pair_distribution_function.x: pair_distribution_function.o mts0_io.o
	${CC} pair_distribution_function.o mts0_io.o -o pair_distribution_function.x

pair_distribution_function.o: pair_distribution_function.cpp
	${CC} -c pair_distribution_function.cpp -Wall -pedantic

mts0_io.o: ../mts0_io/mts0_io.cpp ../mts0_io/mts0_io.h
	${CC} -c ../mts0_io/mts0_io.cpp

clean:	
	rm *.o
