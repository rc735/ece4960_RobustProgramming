COMPILER = g++
COMPILER_FLAGS = -c -g -Wall -Werror

LINKER = g++
LINKER_FLAGS = -g -Wall -Werror

all : main

main : main.cpp matrix_ops.h matrix_ops.cpp param_extraction.h param_extraction.cpp model_transistor.h model_transistor.cpp
	$(COMPILER) $(COMPILER_FLAGS) main.cpp
#matrix_ops.o
	$(COMPILER) $(COMPILER_FLAGS) matrix_ops.cpp
#param_extraction.o
	$(COMPILER) $(COMPILER_FLAGS) param_extraction.cpp
#model_transistors.o
	$(COMPILER) $(COMPILER_FLAGS) model_transistor.cpp
#out
	$(LINKER) main.o matrix_ops.o param_extraction.o model_transistor.o -o out

clean:
	-rm -f *.o out
