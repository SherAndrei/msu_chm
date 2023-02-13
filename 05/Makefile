CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic -Werror=declaration-after-statement -I$(HDR) \
		  -Wmissing-declarations 
LDFLAGS := # -fsanitize=address -fno-omit-frame-pointer -static-libasan
LIBS    := -lm

release: CFLAGS += -O3 -DNDEBUG
release: all

all: Generate.out InterpolationPolynom.out SolveLinearEquations.out

Generate.out: Generate.o ExactSolution.o
	$(CC) $^ -o $@ $(LIBS) $(LDFLAGS) $(CFLAGS)

InterpolationPolynom.out: InterpolationPolynom.o GaussianElimination.o ExactSolution.o
	$(CC) $^ -o $@ $(LIBS) $(LDFLAGS)

SolveLinearEquations.out: SolveLinearEquations.o GaussianElimination.o
	$(CC) $^ -o $@ $(LIBS) $(LDFLAGS) $(CFLAGS)

%.o: %.c Error.h
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean debug release

clean:
	rm -f *.out *.o *.txt

debug: CFLAGS += -g -O0
debug: all
