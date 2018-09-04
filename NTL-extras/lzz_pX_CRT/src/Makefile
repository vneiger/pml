
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
	if [ ! -f ../thresholds_geometric_evaluate_FFT.h ]; then cp -f ../thresholds_geometric_evaluate_FFT.orig $(INC_DIR)/thresholds_geometric_evaluate_FFT.h ; fi
	if [ ! -f ../thresholds_geometric_interpolate_FFT.h ]; then cp -f ../thresholds_geometric_interpolate_FFT.orig $(INC_DIR)/thresholds_geometric_interpolate_FFT.h ; fi
	cp -f ../*.h $(INC_DIR)

clean:
	rm -f $(OBJS)

$(LIB_DIR)/%.o:	%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

