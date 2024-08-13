linalg.o: linalg.cpp linalg.hpp
	g++ -c -O linalg.cpp

test.o: test.cpp linalg.hpp
	g++ -c -O test.cpp

test: linalg.o test.o
	g++ -O -o test linalg.o test.o

clean:
	rm -f *.o test
