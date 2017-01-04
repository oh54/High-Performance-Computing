TARGET	= libmatmult.so
LIBSRCS	= 
LIBOBJS	= ass1_lib.o

OPT	= -g # -O3
PIC	= -fPIC

CXX	= g++
CXXFLAGS= $(OPT) $(PIC) $(XOPTS)

SOFLAGS = -shared 
XLIBS	= 

$(TARGET): $(LIBOBJS)
	$(CXX) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS)

clean:
	@/bin/rm -f core core.* $(LIBOBJS) 

ass1_lib.o : ass1_lib.cpp ass1_lib.h
