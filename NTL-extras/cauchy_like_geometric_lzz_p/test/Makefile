
INC_DIR = ../../include
LIB_DIR = ../../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard test*cpp)
TESTS = $(patsubst %.cpp, u%, $(SRC))
EXEC = $(patsubst %.cpp, %, $(SRC))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lpthread -lgmp -lm
.PHONY: clean

clean:
	rm -f $(EXEC)

all:	$(TESTS)

utest_%:	$(OBJS)
	if test -e test_$*.cpp; then $(CXX) -o test_$* test_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi
