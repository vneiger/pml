
INC_DIR = ../../include
LIB_DIR = ../../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard time*cpp)
EXEC = $(patsubst %.cpp, %, $(SRC))
DATA = $(patsubst %.cpp, %.dat, $(SRC))
RESULTS = $(patsubst %.cpp, %.res, $(SRC))
COMPILE = $(patsubst %.cpp, %.exe, $(patsubst time_%, %, $(SRC)))
RUN = $(patsubst %.cpp, %.dat, $(patsubst time_%, %, $(SRC)))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lpthread -lgmp -lm

.PHONY: clean

clean:
	rm -f $(EXEC) $(DATA) $(RESULTS)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e time_$*.cpp; then $(CXX) -o time_$* time_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi

%.dat: 
	if test -e time_$*; then ./time_$* > time_$*.dat; fi
	if test -e time_$*.awk; then awk -f time_$*.awk time_$*.dat > time_$*.res; fi
	if test -e time_$*.sh; then ./time_$*.sh; fi
