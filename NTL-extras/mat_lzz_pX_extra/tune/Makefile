
INC_DIR = ../../include
LIB_DIR = ../../lib

OBJS = $(wildcard $(LIB_DIR)/*)

SRC = $(wildcard tune*cpp)
EXEC = $(patsubst %.cpp, %, $(SRC))
DATA = $(patsubst %.cpp, %.dat, $(SRC))
RESULTS = $(patsubst %.cpp, %.res, $(SRC))
COMPILE = $(patsubst %.cpp, %.exe, $(patsubst tune_%, %, $(SRC)))
RUN = $(patsubst %.cpp, %.dat, $(patsubst tune_%, %, $(SRC)))

CXX = g++
CXXFLAGS = -Wall -std=c++11 -O3 -march=native -I$(INC_DIR)
LIBS = -lntl -lpthread -lgmp -lm

LBCXXFLAGS = -DHAVE_CONFIG_H -O2 -Wall -g -DNDEBUG -U_LB_DEBUG -DDISABLE_COMMENTATOR -fopenmp -fabi-version=6 -I/usr/local/include -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -mfma
LBLIBS = -lntl -lmpfr -L/usr/local/lib -fopenmp -lopenblas -lpthread -lgfortran -lgivaro -lgmp -lgmpxx -fopenmp

.PHONY: clean

clean:
	rm -f $(EXEC) $(DATA) $(RESULTS)

all:	$(COMPILE)

run:	$(RUN)

%.exe:
	if test -e tune_$*.cpp; then $(CXX) -o tune_$* tune_$*.cpp $(OBJS) $(CXXFLAGS) $(LIBS); fi

%.dat: 
	if test -e tune_$*; then ./tune_$* > tune_$*.dat; fi
	if test -e tune_$*.awk; then \
		if test -e tune_$*.dat; then \
			awk -f tune_$*.awk tune_$*.dat > tune_$*.res;\
		fi; \
	fi
	if test -e tune_$*.sh; then ./tune_$*.sh; fi
