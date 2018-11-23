GFCFLAGS = -O3 -Wall -Wextra


fine.a: fine.f90
	gfortran -fPIC $(GFCFLAGS) -c fine.f90 -I../ffine
	ar rcs fine.a fine.o
	echo "OPEN fine.a\n\
	      ADDLIB ../ffine/fine.a\n\
	      SAVE\n\
	      END" | ar -M
