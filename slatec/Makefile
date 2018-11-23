slatec.a: slatec_src.tgz
	tar -xf slatec_src.tgz
	cd obj; gfortran -std=legacy -fPIC -O3 -c ../src/*.f
	ar rcs slatec.a obj/*.o

slatec_src.tgz:
	wget -N netlib.org/slatec/slatec_src.tgz
