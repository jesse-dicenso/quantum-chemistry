linalg.o: linalg.cpp linalg.hpp
	g++ -O -c linalg.cpp -llapack -lblas

test.o: test.cpp linalg.hpp
	g++ -O -c test.cpp -llapack -lblas

test: linalg.o test.o
	g++ -O -o test linalg.o test.o -llapack -lblas
