CC = mpicc


all: stream_c.exe

stream_c.exe: stream.c
	$(CC) mpibcast_trivial.c -o mpibcast_trivial.exe

clean:
	rm -f mpibcast_trivial.exe *.o

# an example of a more complex build line for the Intel icc compiler
stream.icc: stream.c
	mpicc mpibcast_trivial.c -o mpibcast_trivial
