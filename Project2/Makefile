TARGET = ass2_main

OBJS = ass2_main.o jacobi_seq.o gauss_seidel.o
CXX = suncc
FLAGS = -g -fast -xopenmp

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(FLAGS)

clean:
	@/bin/rm -f *.o

clean_all:
	@/bin/rm -f *.o $(TARGET)


