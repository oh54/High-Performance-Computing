TARGET = ass2_main_jacobi
OBJS = ass2_main.o jacobi_seq.o
CXX = suncc

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) -fast

clean:
	@/bin/rm -f *.o

clean_all:
	@/bin/rm -f *.o $(TARGET)


