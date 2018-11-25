GFCFLAGS = -Wall -Wextra -g
GCCFLAGS = --std=c++17 -Wall -Wextra -g


all: fine.a test

numerics.a: memoise.cpp numerics.f90
	gfortran -fPIC $(GFCFLAGS) -c numerics.f90
	g++ -fPIC $(GCCFLAGS) -c memoise.cpp
	ar rcs numerics.a numerics.o memoise.o

heston.o: numerics.a ../symb_diffs/heston_psi.f90 ../symb_diffs/heston_psi_grad.f90
	gfortran -fPIC $(GFCFLAGS) -I../symb_diffs -c heston.f90

fine.a: heston.o numerics.a ../slatec/slatec.a
	gfortran -fPIC $(GFCFLAGS) -c fine.f90
	ar rcs fine.a fine.o heston.o
	echo "OPEN fine.a\n\
	      ADDLIB numerics.a\n\
	      ADDLIB ../slatec/slatec.a\n\
	      SAVE\n\
	      END" | ar -M

test: fine.a test.f90
	gfortran -o test $(GFCFLAGS) test.f90 fine.a -lstdc++
	./test
