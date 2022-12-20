#Commands
MKDIR := mkdir
RMDIR := rm -rf

#Compiler features
CC      := gcc
CXXFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic \
		  -Wmissing-declarations
LIBS    := -lm

all: Tridiagonal.out Fourier.out

%.out: main.o Task.o %.o
	$(CC) $^ -o $@ $(LIBS)

main.o: main.c Task.h Solve.h
	$(CC) $< -c $(CXXFLAGS)

Task.o: Task.c Task.h
	$(CC) $< -c $(CXXFLAGS)

%.o: %.c Solve.h
	$(CC) $< -c $(CXXFLAGS)

.PHONY: clean debug release

clean:
	$(RMDIR) *.out *.o *.txt *.png

asan: CXXFLAGS += -g -O0 -DNDEBUG -fsanitize=address
asan: LIBS += -static-libasan -fsanitize=address
asan: all

debug: CXXFLAGS += -g -O0
debug: all

release: CXXFLAGS += -O3 -DNDEBUG
release: all
