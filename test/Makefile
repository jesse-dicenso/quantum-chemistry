miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

gf.o: ../libgfs/gf.cpp ../libgfs/gf.hpp ../libmath/miscmath.hpp
	g++ -c -O ../libgfs/gf.cpp

0e.o: ../libint/0e.cpp ../libint/0e.hpp ../libgfs/gf.hpp
	g++ -c -O ../libint/0e.cpp

test.o: test.cpp ../libgfs/gf.hpp ../libint/0e.hpp
	g++ -c -O test.cpp

test: gf.o miscmath.o 0e.o test.o
	g++ -O -o test gf.o miscmath.o 0e.o test.o

clean:
	rm -f *.o
