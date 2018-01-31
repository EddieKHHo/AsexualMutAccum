all: main

main: main.o World.o Popul.o Indiv.o 
	g++ main.o World.o Indiv.o Popul.o -lgsl -lgslcblas -o main
	 
main.o: main.cpp World.h Popul.h Indiv.h IncludeFiles.h 
	g++ -c main.cpp
	
World.o: World.cpp World.h Popul.cpp Popul.h Indiv.h IncludeFiles.h
	g++ -c World.cpp
	 
Popul.o: Popul.cpp Popul.h Indiv.h IncludeFiles.h
	g++ -c Popul.cpp

Indiv.o: Indiv.cpp Indiv.h IncludeFiles.h	 
	g++ -c Indiv.cpp
	
clean:
	 rm main.o main