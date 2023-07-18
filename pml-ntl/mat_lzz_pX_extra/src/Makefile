
INC_DIR = ../../include
LIB_DIR = ../../lib

SOURCES = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp, $(LIB_DIR)/%.o, $(SOURCES))

THRESHOLDS_HEADERS = thresholds_waksman_evaluate thresholds_mp_naive_evaluate thresholds_plain_newton_geometric_mp_inv_trunc thresholds_lmultiplier thresholds_solve_lift

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR) -g
LIBS = -lntl -lgmp -lm
.PHONY: clean

all:	cp $(OBJS)

cp:
	for file in $(THRESHOLDS_HEADERS) ; do \
		if [ ! -f ../$$file.h ]; then cp -f ../$$file.orig $(INC_DIR)/$$file.h ; fi ; \
	done
	cp -f ../*.h $(INC_DIR)

clean:
	rm -f $(OBJS)

reset:	
	for file in $(THRESHOLDS_HEADERS) ; do \
		rm -f ../$$file.h; \
	done


$(LIB_DIR)/%.o:	%.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

