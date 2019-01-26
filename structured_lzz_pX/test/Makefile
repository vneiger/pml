INC_DIR = ../../include
LIB_DIR = ../../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard test*cpp)
EXEC = $(patsubst %.cpp, %, $(SRC))
CHK = $(patsubst %.cpp, %.chk, $(SRC))
COMPILE = $(patsubst %.cpp, %.exe, $(patsubst test_%, %, $(SRC)))
RUN = $(patsubst %.cpp, %.chk, $(patsubst test_%, %, $(SRC)))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lpthread -lgmp -lm

.PHONY: clean

clean:
	rm -f $(EXEC) $(CHK)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e test_$*.cpp; then $(CXX) -o test_$* test_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi

%.chk:
	if test -e test_$*; then ./test_$* > test_$*.chk; fi
