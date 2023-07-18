INC_DIR = ../../include
LIB_DIR = ../../lib
SOURCES = $(wildcard time*c)
OBJS = $(wildcard $(LIB_DIR)/*)

EXEC = $(patsubst %.c, %, $(SOURCES))
DATA = $(patsubst %.c, %.dat, $(SOURCES))
COMPILE = $(patsubst %.c, %.exe, $(patsubst time_%, %, $(SOURCES)))
RUN = $(patsubst %.c, %.dat, $(patsubst time_%, %, $(SOURCES)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lflint -lmpfr -lgmp -lm -lopenblas

.PHONY: clean

clean:
	rm -f $(EXEC) $(DATA)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e time_$*.c; then $(CC) -o time_$* time_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi
