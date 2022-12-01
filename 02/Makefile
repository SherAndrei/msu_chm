#Commands
MKDIR := mkdir
RMDIR := rm -rf

#Compiler features
CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic \
		  -Wmissing-declarations
LIBS    := -lm

all: scheme1.out scheme2.out scheme3.out scheme4.out scheme5.out scheme6.out

scheme%.out: step%.o main.o
	$(CC) $^ -o $@ $(LIBS)

main.o: main.c step.h
	$(CC) $< -c $(CFLAGS)

step%.o: step%.c step.h
	$(CC) $< -c $(CFLAGS)

.PHONY: clean debug release

clean:
	$(RMDIR) *.gch *.o *.out

asan: CFLAGS += -g -O0 -DNDEBUG -fsanitize=address -static-libasan
asan: all

debug: CFLAGS += -g -O0
debug: all

release: CFLAGS += -O3 -DNDEBUG
release: all
