#CC = icpc -I../mts0_io -O3
CC = g++ -std=c++11
PROG = measure_density_voro
VORO = -I/usr/local/include/voro++ -L/usr/local/lib -lvoro++
MTS0 = -I../mts0_io

${PROG}.x: ${PROG}.o mts0_io.o
	${CC} ${PROG}.o mts0_io.o ${VORO} -o ${PROG}.x

${PROG}.o: ${PROG}.cpp
	${CC} -c ${PROG}.cpp ${VORO} ${MTS0} -Wall -pedantic

mts0_io.o: ../mts0_io/mts0_io.cpp ../mts0_io/mts0_io.h
	${CC} -c ../mts0_io/mts0_io.cpp ${MTS0}

clean:	
	rm *.o
