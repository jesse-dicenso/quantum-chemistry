miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

gf.o: ../libgf/gf.cpp ../libgf/gf.hpp ../libmath/miscmath.hpp
	g++ -c -O ../libgf/gf.cpp

1e.o: ../libint/1e.cpp ../libint/1e.hpp ../libgf/gf.hpp
	g++ -c -O ../libint/1e.cpp

2e.o: ../libint/2e.cpp ../libint/2e.hpp ../libint/1e.hpp
	g++ -c -O ../libint/2e.cpp

test.o: test.cpp ../libint/2e.hpp 
	g++ -c -O test.cpp

test: miscmath.o gf.o 1e.o 2e.o test.o
	g++ -O -o test miscmath.o gf.o 1e.o 2e.o test.o

clean:
	rm -f *.o test
