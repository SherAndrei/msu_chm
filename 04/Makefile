#Commands
MKDIR := mkdir
RMDIR := rm -rf

#Compiler features
CC      := g++
CXXFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic \
		  -Wmissing-declarations -std=c++14
LIBS    := -lm

all: Tridiagonal.out Fourier.out

%.out: main.o Task.o Vector.o %.o
	$(CC) $^ -o $@ $(LIBS)

main.o: main.cpp Task.h Solve.h Vector.h
	$(CC) $< -c $(CXXFLAGS)

Vector.o: Vector.cpp Vector.h
	$(CC) $< -c $(CXXFLAGS)

Task.o: Task.cpp Task.h
	$(CC) $< -c $(CXXFLAGS)

%.o: %.cpp Solve.h Vector.h
	$(CC) $< -c $(CXXFLAGS)

.PHONY: clean debug release

clean:
	$(RMDIR) *.out *.o *.txt

asan: CXXFLAGS += -g -O0 -DNDEBUG -fsanitize=address
asan: LIBS += -static-libasan -fsanitize=address
asan: all

debug: CXXFLAGS += -g -O0
debug: all

release: CXXFLAGS += -O3 -DNDEBUG
release: all
