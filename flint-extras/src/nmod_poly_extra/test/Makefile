SYSINC_DIR = ../../include
INC_DIR = ../../include
SYSLIB_DIR = /usr/local/lib
LIB_DIR = ../../lib
OBJS = $(wildcard $(LIB_DIR)/*.o)
SRC = $(wildcard test*c)
EXEC = $(patsubst %.c, %, $(SRC))
CHK = $(patsubst %.c, %.chk, $(patsubst test_%, %, $(SRC)))
COMPILE = $(patsubst %.c, %.exe, $(patsubst test_%, %, $(SRC)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I$(SYSINC_DIR) -I/usr/local/include -I/opt/homebrew/include -L$(LIB_DIR) -L$(SYSLIB_DIR) -L/opt/homebrew/lib
LIBS = -lflint -lmpfr -lgmp -lm -lpthread

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(CHK)

%.exe:
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi

%.chk:
	if test -e test_$*; then ./test_$* > $*.chk; fi
