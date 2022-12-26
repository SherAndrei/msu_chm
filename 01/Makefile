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

0011.out: main.c
	$(CC) $^ -o $@ $(CFLAGS) $(LIBS)

.PHONY: clean debug release

clean:
	$(RMDIR) *.out *.o

asan: CFLAGS += -g -O0 -DNDEBUG -fsanitize=address -static-libasan
asan: all

debug: CFLAGS += -g -O0
debug: all

release: CFLAGS += -O3 -DNDEBUG
release: all
