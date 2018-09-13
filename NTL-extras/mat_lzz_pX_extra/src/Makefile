
INC_DIR = ../../include
LIB_DIR = ../../lib

SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(LIB_DIR)/%.o, $(SOURCES))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR) -g
LIBS = -lntl -lgmp -lm
.PHONY: clean

all:	cp $(OBJS)

cp:
	if [ ! -f ../thresholds_waksman_evaluate.h ]; then cp -f ../thresholds_waksman_evaluate.orig $(INC_DIR)/thresholds_waksman_evaluate.h ; fi
	if [ ! -f ../thresholds_mp_naive_evaluate.h ]; then cp -f ../thresholds_mp_naive_evaluate.orig $(INC_DIR)/thresholds_mp_naive_evaluate.h ; fi

	cp -f ../*.h $(INC_DIR)

clean:
	rm -f $(OBJS)

$(LIB_DIR)/%.o:	%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

