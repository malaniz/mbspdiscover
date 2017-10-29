MCBSPLIB=../libs/mcbsp/lib/libmcbsp1.2.0.a
MCBSPINCLUDE=../libs/mcbsp/
HWLOCLIB=../libs/hwloc/libhwloc.so
HWLOCINCLUDE=../libs/hwloc/include



#CC=icc -mmic
CC=gcc
CFLAGS=-I $(MCBSPINCLUDE) -I $(HWLOCINCLUDE)
LFLAGS=$(MCBSPLIB) -pthread -lm -lrt $(HWLOCLIB)
OBJ= mbspbench.o mbsp-discover.o mbsputil.o

all: mbspbench

mbspbench: $(OBJ)
	$(CC) $(CFLAGS) -o mbspbench $(OBJ) $(LFLAGS)

mbspbench.o: mbspbench.c
	$(CC) $(CFLAGS) -c mbspbench.c

mbsp-discover.o: mbsp-discover.c mbsp-discover.h
	$(CC) $(CFLAGS) -c mbsp-discover.c

mbsputil.o: mbsputil.c  mbsputil.h
	$(CC) $(CFLAGS) -c mbsputil.c

clean:
	rm -f *.o mbspbench
