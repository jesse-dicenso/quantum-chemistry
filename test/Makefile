miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -Wall -c -O ../libmath/miscmath.cpp

linalg.o: ../libmath/linalg.cpp ../libmath/linalg.hpp
	g++ -Wall -c -O ../libmath/linalg.cpp

gf.o: ../libgf/gf.cpp ../libgf/gf.hpp ../libmath/miscmath.hpp
	g++ -Wall -c -O ../libgf/gf.cpp

1e.o: ../libint/1e.cpp ../libint/1e.hpp ../libgf/gf.hpp ../libmath/linalg.hpp
	g++ -Wall -c -O ../libint/1e.cpp

2e.o: ../libint/2e.cpp ../libint/2e.hpp ../libint/1e.hpp
	g++ -Wall -c -O ../libint/2e.cpp

mol.o: ../libmol/mol.cpp ../libmol/mol.hpp ../libgf/gf.hpp
	g++ -Wall -c -O ../libmol/mol.cpp

test.o: test.cpp ../libint/2e.hpp ../libmol/mol.hpp
	g++ -Wall -c -O test.cpp

test: miscmath.o linalg.o gf.o 1e.o 2e.o mol.o test.o
	g++ -Wall -O -o test miscmath.o linalg.o gf.o 1e.o 2e.o mol.o test.o

clean:
	rm -f *.o
