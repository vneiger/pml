
SYSINC_DIR = /usr/local/include
INC_DIR = ../../include
INC_MAC = /opt/homebrew/include

SYSLIB_DIR = /usr/local/lib
LIB_DIR = ../../lib
LIB_MAC = /opt/homebrew/lib

OBJS = $(wildcard $(LIB_DIR)/*.o)
SRC = $(wildcard test*c)
EXEC = $(patsubst %.c, %, $(SRC))
CHK = $(patsubst %.c, %.chk, $(patsubst test_%, %, $(SRC)))
COMPILE = $(patsubst %.c, %.exe, $(patsubst test_%, %, $(SRC)))

CC = gcc
CFLAGS = -Wall -O3 -march=native -I$(INC_DIR) -I$(INC_MAC) -I$(SYSINC_DIR) -L$(LIB_DIR) -L$(LIB_MAC) -L$(SYSLIB_DIR) 
LIBS = -lflint -lmpfr -lgmp -lm -lpthread -lpml

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(CHK)

%.exe:
	if test -e test_$*.c; then $(CC) -o test_$* test_$*.c $(OBJS) $(CFLAGS) $(LIBS); fi

%.chk:
	if test -e test_$*; then ./test_$* > $*.chk; fi
