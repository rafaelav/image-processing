SHELL = /bin/sh
CC    = mpicc
MAKE  = make
LIBS  = -lm -lmpi
MPIRUN= mpirun
NP_ARG= -np 

all: prelucrare_imagine

SOURCES = prelucrare_imagine.c 
PROGRAMS = prelucrare_imagine

clean: 
	rm -f *.o ${PROGRAMS} *~ PI*

prelucrare_imagine: ${SOURCES}
	${CC} ${SOURCES} ${LIBS} -o exec

run: prelucrare_imagine 
	${MPIRUN} ${NP_ARG} 4 ./exec 2 picture.pgm 1 1 1 outpicture.pgm

sources: @echo ${SOURCES} 










