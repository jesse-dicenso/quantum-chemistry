gf.o : gf.cpp gf.hpp
	g++ -c -O gf.cpp

test.o : test.cpp gf.hpp
	g++ -c -O test.cpp

test : gf.o test.o
	g++ -O -o test gf.o test.o
