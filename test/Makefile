miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

gf.o: ../libgf/gf.cpp ../libgf/gf.hpp ../libmath/miscmath.hpp
	g++ -c -O ../libgf/gf.cpp

1e.o: ../libint/1e.cpp ../libint/1e.hpp ../libgf/gf.hpp
	g++ -c -O ../libint/1e.cpp

test.o: test.cpp ../libint/1e.hpp 
	g++ -c -O test.cpp

test: gf.o miscmath.o 1e.o test.o
	g++ -O -o test gf.o miscmath.o 1e.o test.o

clean:
	rm -f *.o test
