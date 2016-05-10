COMPILER = g++
COMPILER_FLAGS = -c -g -Wall
INC = -I/program4 -I/program3

LINKER = g++

all : main

main : main.cpp
	cp define.h program3/defines.h
	cp define.h program4/define.h
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program3/pg3_model_transistor.cpp
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program3/pg3_param_extraction.cpp
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program3/pg3_matrix_ops.cpp
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program4/pg4_methods.cpp
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program4/pg4_solver.cpp
	$(COMPILER) $(COMPILER_FLAGS) $(INC) program4/pg4_param_extraction.cpp

	$(COMPILER) $(COMPILER_FLAGS) main.cpp
	$(LINKER) $(INC) *.o -o out

clean :
	rm -rf *.o out
