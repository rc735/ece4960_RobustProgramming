COMPILER = g++
COMPILER_FLAGS = -c -g -Wall

LINKER = g++

all : main

main : main.cpp methods.cpp solver.cpp param_extraction.cpp matrix_ops.cpp
	$(COMPILER) $(COMPILER_FLAGS) main.cpp
	$(COMPILER) $(COMPILER_FLAGS) methods.cpp
	$(COMPILER) $(COMPILER_FLAGS) solver.cpp
	$(COMPILER) $(COMPILER_FLAGS) param_extraction.cpp
	$(COMPILER) $(COMPILER_FLAGS) matrix_ops.cpp
	$(LINKER) *.o -o out

matrix_ops : matrix_ops_test.cpp matrix_ops.cpp
	$(COMPILER) matrix_ops_test.cpp matrix_ops.cpp -o matout

clean : 
	rm -rf *.o out matout
