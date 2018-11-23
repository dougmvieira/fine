GCCFLAGS = -O3 -fPIC -shared --std=c++17 -Wall -Wextra
PYINCLUDE = $(shell python3 -c "import sysconfig; print(sysconfig.get_path('include', 'posix_prefix'))")
INCLUDES = -I/usr/include/eigen3 -I$(PYINCLUDE) -I../cppfine


fine.so: fine.cpp ../cppfine/fine.hpp
	g++ $(GCCFLAGS) fine.cpp ../cppfine/fine.a -o fine.so $(INCLUDES) -lgfortran
	python3 test.py
