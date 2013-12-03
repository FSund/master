CC = g++ -I../mts0_io -O3
PROG = tag_atoms

${PROG}.x: ${PROG}.o mts0_io.o
	${CC} ${PROG}.o mts0_io.o -o ${PROG}.x

${PROG}.o: ${PROG}.cpp
	${CC} -c ${PROG}.cpp -Wall -pedantic

mts0_io.o: ../mts0_io/mts0_io.cpp ../mts0_io/mts0_io.h
	${CC} -c ../mts0_io/mts0_io.cpp

clean:  
	rm *.o
