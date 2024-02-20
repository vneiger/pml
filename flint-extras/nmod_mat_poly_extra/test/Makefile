INC_DIR = ../../include
LIB_DIR = ../../lib
OBJS = $(wildcard $(LIB_DIR)/*)
SRC = $(wildcard test*c)
EXEC = $(patsubst %.c, %, $(SRC))
CHK = $(patsubst %.c, %.chk, $(SRC))
COMPILE = $(patsubst %.c, %.exe, $(patsubst test_%, %, $(SRC)))
RUN = $(patsubst %.c, %.chk, $(patsubst test_%, %, $(SRC)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR)
LIBS = -lflint -lmpfr -lgmp -lm

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi

%.chk:
	if test -e test_$*; then ./test_$* > test_$*.chk; fi
