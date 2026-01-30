CXX := g++
CXXFLAGS := -O2 -Wno-c++17-extensions #-fopenmp

SRC_DIRS := . ../libaux ../libdft ../libgf ../libgrid ../libint ../libmath ../libmol ../libscf
SRC := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))
OBJ := $(SRC:.cpp=.o)
TARGET := QC-EXEC

LIBS := -lblas -llapack #-lgomp

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

.PHONY: all clean
