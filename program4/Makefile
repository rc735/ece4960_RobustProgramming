COMPILER = g++
COMPILER_FLAGS = -c -g -Wall

LINKER = g++

all : main

main : main.cpp methods.cpp solver.cpp
	$(COMPILER) $(COMPILER_FLAGS) main.cpp
	$(COMPILER) $(COMPILER_FLAGS) methods.cpp
	$(COMPILER) $(COMPILER_FLAGS) solver.cpp
	$(LINKER) main.o methods.o solver.o -o out

clean : 
	rm -rf *.o out
