CXX = g++
CXXFLAGS = -g -O -fsanitize=address
LDFLAGS = -fsanitize=address -llapack -lblas

SOURCES = ../libmath/miscmath.cpp ../libmath/linalg.cpp ../libgf/gf.cpp ../libint/1e.cpp ../libint/2e.cpp ../libmol/mol.cpp hartree-fock.cpp
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = hartree-fock

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS)

%.o: %.cpp
	$(CXX) $(LDFLAGS) -c $< -o $@

clean:
	rm -f *.o
