TARGET = ass2_main
OBJS = ass2_main.o jacobi_seq.o gauss_seidel.o
CC = suncc
CFLAGS = -fast -xopenmp


$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) 

clean:
	@/bin/rm -f *.o

clean_all:
	@/bin/rm -f *.o $(TARGET)

ass2_main.o : ass2_main.c ass2_lib.h
jacobi_seq.o : jacobi_seq.c
gauss_seidel.o : gauss_seidel.c

