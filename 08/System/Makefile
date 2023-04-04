CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic -Werror=declaration-after-statement \
		  -Wmissing-declarations -Werror=strict-prototypes
LDFLAGS := -lm

release: CFLAGS += -O3 -DNDEBUG
release: all

all: Newton.out InvertMatrix.out DiffEquation.out

InvertMatrix.out: InvertMatrix.o GaussianInversion.o
	$(CC) $^ -o $@ $(LDFLAGS)

Newton.out: Newton.o GaussianInversion.o Root.o
	$(CC) $^ -o $@ $(LDFLAGS)

DiffEquation.out:  DiffEquation.o Root.o GaussianInversion.o
	$(CC) $^ -o $@ $(LDFLAGS)

Root.o: Root.c Root.h GaussianInversion.h Error.h
	$(CC) -c $< -o $@ $(CFLAGS)

%.o: %.c Error.h Root.h
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean debug release asan

clean:
	rm -f *.out *.o *.txt

debug: CFLAGS += -g -O0
debug: all

asan: CFLAGS += -g -O0 -fsanitize=address
asan: LDFLAGS += -fsanitize=address
asan: all
