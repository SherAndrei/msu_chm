#Compiler features
CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic \
		  -Wmissing-declarations
LIBS    := -lm

release: CFLAGS += -O3 -DNDEBUG
release: all

asan: CFLAGS += -fsanitize=address -g -O0
asan: LIBS += -fsanitize=address
asan: all

all: \
	AlphaStable.out \
	Alternating.out \
	LebedevMethod.out \
	ExplicitEuler.out \
	ImplicitAdams.out \
	ImplicitEuler.out \
	NotAlphaStable.out \
	SeeminglyStable.out \
	TestScheme.out \

LebedevMethod.out: main.o LebedevMethod.o
	$(CC) $^ -o $@ $(LIBS)

TestScheme.out: main.o TestScheme.o SingleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

NotAlphaStable.out: ExactSolution.o main.o NotAlphaStable.o DoubleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

SeeminglyStable.out: ExactSolution.o main.o SeeminglyStable.o DoubleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

AlphaStable.out: ExactSolution.o main.o AlphaStable.o DoubleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

Alternating.out: ExactSolution.o main.o Alternating.o DoubleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

%.out: ExactSolution.o main.o %.o SingleStepScheme.o
	$(CC) $^ -o $@ $(LIBS)

%.o: %.c ExactSolution.h
	$(CC) $< -c $(CFLAGS)

.PHONY: clean debug release

clean:
	rm -rf *.gch *.o *.out

debug: CFLAGS += -g -O0
debug: all
