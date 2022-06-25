INC_DIR = ../../include
LIB_DIR = ../../lib
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard ../*.h)
OBJS = $(patsubst %.c, $(LIB_DIR)/%.o, $(SOURCES))
INSTALLED_HEADERS = $(patsubst ../%.h, $(INC_DIR)/%.h, $(HEADERS))

CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native -I$(INC_DIR)

.PHONY: clean

all:    cp $(OBJS)

cp:
	cp -f $(HEADERS) $(INC_DIR)

clean:
	rm -f $(INSTALLED_HEADERS)
	rm -f $(OBJS)

$(LIB_DIR)/%.o:       %.c
	$(CC) -c $< $(CFLAGS) -o $@
