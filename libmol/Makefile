mol.o: mol.cpp mol.hpp
	g++ -c mol.cpp

test.o: test.cpp mol.hpp
	g++ -c test.cpp

test: mol.o test.o
	g++ -O -o test mol.o test.o

clean:
	rm -f *.o
