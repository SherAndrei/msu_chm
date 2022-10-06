#Commands
MKDIR := mkdir
RMDIR := rm -rf

#Compiler features
CC      := ../bgch_cmpl/g++
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic \
		  -Wmissing-declarations -std=c++14
LIBS    := -lm

NAME := 0011
EXE := ./$(NAME).out

$(EXE): $(NAME).cpp
	$(CC) $^ -o $@ $(CFLAGS) $(LIBS)

.PHONY: clean debug release

clean:
	$(RMDIR) $(EXE)

asan: CFLAGS += -g -O0 -DNDEBUG -fsanitize=address -static-libasan
asan: $(EXE)

debug: CFLAGS += -g -O0
debug: $(EXE)

release: CFLAGS += -O3 -DNDEBUG
release: $(EXE)
