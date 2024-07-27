miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

gf.o: gf.cpp gf.hpp ../libmath/miscmath.hpp
	g++ -c -O gf.cpp

test.o: test.cpp gf.hpp ../libmath/miscmath.hpp
	g++ -c -O test.cpp

test: gf.o miscmath.o test.o
	g++ -O -o test gf.o miscmath.o test.o

clean:
	rm -f *.o
