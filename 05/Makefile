CC      := gcc
CFLAGS  := -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align \
		  -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security \
		  -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long  \
		  -Wcast-qual -Wno-suggest-attribute=format -Wpedantic -Werror=declaration-after-statement -I$(HDR) \
		  -Wmissing-declarations 
LDFLAGS := # -fsanitize=address -fno-omit-frame-pointer -static-libasan
LIBS    := -lm

all: Generate.out InterpolationPolynom.out

Generate.out: Generate.c
	$(CC) $^ -o $@ $(LIBS) $(LDFLAGS) $(CFLAGS)

InterpolationPolynom.out: InterpolationPolynom.o Matrix.o
	$(CC) $^ -o $@ $(LIBS) $(LDFLAGS)

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean debug release

clean:
	rm -f *.out *.o
	ls *.txt 2>/dev/null && rm -i *.txt

debug: CFLAGS += -g -O0
debug: all

release: CFLAGS += -O3 -DNDEBUG
release: all
