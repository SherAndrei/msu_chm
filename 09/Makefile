CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic -Werror=declaration-after-statement \
		  -Wmissing-declarations -Werror=strict-prototypes
LDFLAGS := -lm

release: CFLAGS += -O3 -DNDEBUG
release: all

all: Generate.out RBFInterpolator.out

Generate.out: Generate.o ExactSolution.o
	$(CC) $^ -o $@ $(LDFLAGS)

RBFInterpolator.out: RBFInterpolator.o GaussianElimination.o KNearestNeighbors.o
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.c Error.h
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean debug release asan

clean:
	rm -f *.out *.o *.txt *.gch *.png

debug: CFLAGS += -g -O0
debug: all

asan: CFLAGS += -g -O0 -fsanitize=address
asan: LDFLAGS += -fsanitize=address
asan: all
