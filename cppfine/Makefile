GCCFLAGS = --std=c++17 -O3 -Wall -Wextra


fine.a: test
	ar rcs fine.a fine.o
	echo "OPEN fine.a\n\
	      ADDLIB ../cfine/fine.a\n\
	      SAVE\n\
	      END" | ar -M

test: fine.o ../cfine/fine.a
	g++ $(GCCFLAGS) -o test test.cpp fine.o ../cfine/fine.a -lgfortran
	./test

fine.o: fine.cpp fine.hpp
	g++ -fPIC $(GCCFLAGS) -c fine.cpp
