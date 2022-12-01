#Commands
MKDIR := mkdir
RMDIR := rm -rf

#Compiler features
CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic -Werror=declaration-after-statement -I$(HDR) \
		  -Wmissing-declarations
LIBS    := -lm

all: exp.out sin_and_cos.out cubic.out new.out

%.out: step.o main.o %.o
	$(CC) $^ -o $@ $(LIBS)

main.o: main.c step.h system_of_equations.h
	$(CC) $< -c $(CFLAGS)

new.o: new.c system_of_equations.h
	$(CC) $< -c $(CFLAGS)

exp.o: exp.c system_of_equations.h
	$(CC) $< -c $(CFLAGS)

sin_and_cos.o: sin_and_cos.c system_of_equations.h
	$(CC) $< -c $(CFLAGS)

cubic.o: cubic.c system_of_equations.h
	$(CC) $< -c $(CFLAGS)

step.o: step.c step.h
	$(CC) $< -c $(CFLAGS)

.PHONY: clean debug release

clean:
	$(RMDIR) *.gch *.o *.out

asan: CFLAGS += -g -O0 -DNDEBUG -fsanitize=address
asan: all

debug: CFLAGS += -g -O0
debug: all

release: CFLAGS += -O3 -DNDEBUG
release: all
