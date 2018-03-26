#
# Compile X-Score on a Unix/Linux machine with the GNU C++ compiler 
#
CC = g++                           # apply the GNU C++ compiler
#
OBJECT = main.o \
         basic.o \
         score.o \
         parameter.o \
         input.o \
         ligand.o \
         logp.o \
         misc.o \
         molecule.o \
         population.o \
         protein.o \
         utilities.o 
#
LIB = -lm
#	 
xscore: $(OBJECT)
	$(CC) $(OBJECT) $(LIB) -o xscore 
#
main.o: main.c xtool.h
basic.o: basic.c xtool.h
score.o: score.c xtool.h
parameter.o: parameter.c xtool.h
input.o: input.c xtool.h
ligand.o: ligand.c xtool.h
logp.o: logp.c xtool.h
misc.o: misc.c xtool.h
molecule.o: molecule.c xtool.h
population.o: population.c xtool.h
protein.o: protein.c xtool.h
utilities.o: utilities.c xtool.h
#
