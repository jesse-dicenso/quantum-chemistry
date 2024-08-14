miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

linalg.o: ../libmath/linalg.cpp ../libmath/linalg.hpp
	g++ -c -O ../libmath/linalg.cpp

gf.o: ../libgf/gf.cpp ../libgf/gf.hpp ../libmath/miscmath.hpp
	g++ -c -O ../libgf/gf.cpp

1e.o: ../libint/1e.cpp ../libint/1e.hpp ../libgf/gf.hpp
	g++ -c -O ../libint/1e.cpp

2e.o: ../libint/2e.cpp ../libint/2e.hpp ../libint/1e.hpp
	g++ -c -O ../libint/2e.cpp

buildmatrix.o: ../libint/buildmatrix.cpp ../libint/buildmatrix.hpp
	g++ -c -O ../libint/buildmatrix.cpp

test.o: test.cpp ../libint/2e.hpp 
	g++ -c -O test.cpp

test: miscmath.o linalg.o gf.o 1e.o 2e.o buildmatrix.o test.o
	g++ -O -o test miscmath.o linalg.o gf.o 1e.o 2e.o buildmatrix.o test.o

clean:
	rm -f *.o
