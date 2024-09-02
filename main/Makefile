miscmath.o: ../libmath/miscmath.cpp ../libmath/miscmath.hpp
	g++ -c -O ../libmath/miscmath.cpp

linalg.o: ../libmath/linalg.cpp ../libmath/linalg.hpp
	g++ -c -O ../libmath/linalg.cpp -llapack -lblas

gf.o: ../libgf/gf.cpp ../libgf/gf.hpp ../libmath/miscmath.hpp
	g++ -c -O ../libgf/gf.cpp

1e.o: ../libint/1e.cpp ../libint/1e.hpp ../libgf/gf.hpp ../libmath/linalg.hpp
	g++ -c -O ../libint/1e.cpp

2e.o: ../libint/2e.cpp ../libint/2e.hpp ../libint/1e.hpp
	g++ -c -O ../libint/2e.cpp

mol.o: ../libmol/mol.cpp ../libmol/mol.hpp ../libgf/gf.hpp
	g++ -c -O ../libmol/mol.cpp

scfgen.o: ../libscf/scfgen.cpp ../libscf/scfgen.hpp ../libint/2e.hpp ../libmath/linalg.hpp ../libmol/mol.hpp
	g++ -c -O ../libscf/scfgen.cpp

scfrhf.o: ../libscf/scfrhf.cpp ../libscf/scfrhf.hpp ../libint/2e.hpp ../libmath/linalg.hpp ../libmol/mol.hpp
	g++ -c -O ../libscf/scfrhf.cpp

scfuhf.o: ../libscf/scfuhf.cpp ../libscf/scfuhf.hpp ../libint/2e.hpp ../libmath/linalg.hpp ../libmol/mol.hpp
	g++ -c -O ../libscf/scfuhf.cpp

hartree-fock.o: hartree-fock.cpp hartree-fock.hpp ../libscf/scfgen.hpp ../libscf/scfrhf.hpp ../libscf/scfuhf.hpp
	g++ -c -O hartree-fock.cpp

hartree-fock: miscmath.o linalg.o gf.o 1e.o 2e.o mol.o scfgen.o scfrhf.o scfuhf.o hartree-fock.o
	g++ -O -o hartree-fock miscmath.o linalg.o gf.o 1e.o 2e.o mol.o scfgen.o scfrhf.o scfuhf.o hartree-fock.o -llapack -lblas

clean:
	rm -f *.o
