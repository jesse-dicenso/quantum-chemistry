CXX := g++
CXXFLAGS := -O

SRC_DIRS := . ../libdft ../libgf ../libgrid ../libint ../libmath ../libmol ../libscf
SRC := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))
OBJ := $(patsubst %.cpp, build/%.o, $(SRC))
TARGET := QC-EXEC

LIBS := -lblas -llapack

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LIBS) -o $@

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

.PHONY: all clean
