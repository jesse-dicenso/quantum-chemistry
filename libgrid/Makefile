test.o: test.cpp quadrature.hpp
	g++ -c -O test.cpp

quadrature.o: quadrature.cpp quadrature.hpp
	g++ -c -O quadrature.cpp

test: test.o quadrature.o
	g++ -O -o test test.o quadrature.o

clean:
	rm -f *.o
