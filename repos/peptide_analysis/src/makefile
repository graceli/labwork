# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall -g
#-g allo executable to be used with GDB 

all: peptide_analysis

peptide_analysis: peptide_analysis_all.o Inositol.o Peptide.o PepGroup.o HB.o Water.o
	$(CC) peptide_analysis_all.o Inositol.o Peptide.o PepGroup.o HB.o Water.o -o ../exe/peptide_analysis_all

PepGroup.o: PepGroup.cpp PepGroup.h
	$(CC) $(CFLAGS) PepGroup.cpp

peptide_analysis_all.o: peptide_analysis_all.cpp FileRead.h
	$(CC) $(CFLAGS) peptide_analysis_all.cpp

Inositol.o: Inositol.cpp Inositol.h HB.h
	$(CC) $(CFLAGS) Inositol.cpp

Peptide.o: Peptide.cpp Peptide.h
	$(CC) $(CFLAGS) Peptide.cpp

HB.o: HB.cpp HB.h
	$(CC) $(CFLAGS) HB.cpp

Water.o: Water.cpp Water.h
	$(CC) $(CFLAGS) Water.cpp
	
clean:
	rm *.o ../exe/peptide_analysis_all
