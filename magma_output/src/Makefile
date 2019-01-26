
INC_DIR = ../../include
LIB_DIR = ../../lib

SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(LIB_DIR)/%.o, $(SOURCES))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lgmp -lm
.PHONY: clean

all:	cp $(OBJS)

cp:
	cp -f ../*.h $(INC_DIR)

clean:
	rm -f $(OBJS)

$(LIB_DIR)/%.o:	%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

