CC = g++ -I../mts0_io -O3
PROG = get_xyz

${PROG}.x: ${PROG}.o mts0_io.o
	${CC} ${PROG}.o mts0_io.o -o ${PROG}.x

${PROG}.o: ${PROG}.cpp
	${CC} -c ${PROG}.cpp

mts0_io.o: ../mts0_io/mts0_io.cpp ../mts0_io/mts0_io.h
	${CC} -c ../mts0_io/mts0_io.cpp

clean:
	rm *.o
