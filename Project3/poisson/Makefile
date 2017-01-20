TARGET	 = poisson
SRCS	 = 
OBJS	 = ass2_main.o jacobis.o uo_updates.o

OPT      = -g
PIC 	 =
XOPTS	 = -Xptxas=-v
ARCH 	 = -arch=sm_30
OMP  	 = -fopenmp

CXX	 = nvcc
CXXFLAGS = --compiler-options "$(OPT) $(PIC) $(OMP)" $(ARCH) $(XOPTS) 

CUDA_PATH ?= /appl/cuda/8.0
INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/common/inc

SOFLAGS =
XLIBS	= 

$(TARGET): $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS) 

.SUFFIXES: .cu
.cu.o:
	$(CXX) -o $*.o -c $*.cu $(CXXFLAGS) $(SOFLAGS) $(INCLUDES)

clean:
	@/bin/rm -f $(TARGET) $(OBJS) 


ass2_main.o : ass2_main.cu ass2_lib.h
jacobis.o : jacobis.cu
uo_updates.o : uo_updates.cu
