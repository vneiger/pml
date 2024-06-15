SYSINC_DIR = /usr/local/include
INC_DIR = ../../include
SYSLIB_DIR = /usr/local/lib
LIB_DIR = ../../lib
SOURCES = $(wildcard time*c)
OBJS = $(wildcard $(LIB_DIR)/*)
EXEC = $(patsubst %.c, %, $(SOURCES))
DATA = $(patsubst %.c, %.dat, $(SOURCES))
COMPILE = $(patsubst %.c, %.exe, $(patsubst time_%, %, $(SOURCES)))
RUN = $(patsubst %.c, %.dat, $(patsubst time_%, %, $(SOURCES)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I$(SYSINC_DIR) -L$(LIB_DIR) -L$(SYSLIB_DIR) 
#CFLAGS = -fverbose-asm -S -Wall -O3 -march=native -I$(INC_DIR) -I$(SYSINC_DIR) -L$(LIB_DIR) -L$(SYSLIB_DIR) 
LIBS = -lflint -lmpfr -lgmp -lm -lpthread

.PHONY: clean

dot-clang:
	clang -c -o dot.o dot.c -Wall -Wextra -O3 -march=native
	clang -c -o time_dot.o time_dot.c -Wall -Wextra -O3 -march=native
	clang -o time_dot time_dot.o dot.o -lflint

dot-gcc:
	gcc -c -o dot.o dot.c -Wall -Wextra -O3 -march=native
	gcc -c -o time_dot.o time_dot.c -Wall -Wextra -O3 -march=native
	gcc -o time_dot time_dot.o dot.o -lflint

clean:
	rm -f $(EXEC) $(DATA)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e time_$*.c; then $(CC) -o time_$* time_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi
