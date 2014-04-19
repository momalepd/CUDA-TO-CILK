CXX=icpc
BITS=-m64
INCLUDES=-I /home/prashant/momalepd/Codes/CUDA-to-Cilk-Plus/include/
CXXFLAGS=-O2 $(INCLUDES) $(BITS)

# If you are using a version of binutils prior to 2.17, 
# uncomment following line to work around a known issue. 
# please reference the Release Notes for more information. 
#LIBS=-lcilkrts

PROGS= _main

all:: $(PROGS)

test: test.o
	$(CXX) $(BITS) -tbb -o $@ _main.o $(LIBS)


clean::
	-rm -f $(PROGS) *.o
