#CC = icpc -I../mts0_io -O3
CC = g++ -std=c++11
PROG = measure_density_water
VORO = -I/usr/local/include/voro++ -L/usr/local/lib -lvoro++
MTS0 = -I../mts0_io
AVG = measure_average_density_of_atom_types
AVG_WATER = measure_average_density_of_water

all: $(PROG).x $(AVG).x $(AVG_WATER).x

average: $(AVG).x

average_water: $(AVG_WATER).x

${PROG}.x: ${PROG}.o mts0_io.o
	${CC} ${PROG}.o mts0_io.o ${VORO} -o ${PROG}.x

${PROG}.o: ${PROG}.cpp
	${CC} -c ${PROG}.cpp ${VORO} ${MTS0} -Wall -pedantic

${AVG}.x: ${AVG}.o mts0_io.o
	${CC} ${AVG}.o mts0_io.o ${VORO} -o ${AVG}.x

${AVG}.o: ${AVG}.cpp
	${CC} -c ${AVG}.cpp ${VORO} ${MTS0} -Wall -pedantic

${AVG_WATER}.x: ${AVG_WATER}.o mts0_io.o
	${CC} ${AVG_WATER}.o mts0_io.o ${VORO} -o ${AVG_WATER}.x

${AVG_WATER}.o: ${AVG_WATER}.cpp
	${CC} -c ${AVG_WATER}.cpp ${VORO} ${MTS0} -Wall -pedantic

mts0_io.o: ../mts0_io/mts0_io.cpp ../mts0_io/mts0_io.h
	${CC} -c ../mts0_io/mts0_io.cpp ${MTS0}

clean:	
	rm *.o
