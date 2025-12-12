CXX := clang++
CXXFLAGS := -O2

SRC_DIRS := . ../libaux ../libdft ../libgf ../libgrid ../libint ../libmath ../libmol ../libscf
SRC := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))
OBJ := $(SRC:.cpp=.o)
TARGET := QC-EXEC

LIBS := -lblas -llapack

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

.PHONY: all clean
