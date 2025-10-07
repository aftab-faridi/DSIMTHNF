FC = gfortran
CC = gcc

FFLAGS = -O2 -c
CFLAGS = -O2 -c
LDFLAGS =

OBJS = main.o dsimthnf.o

all: run

main.o: main.c
	$(CC) $(CFLAGS) main.c

dsimthnf.o: dsimthnf.f
	$(FC) $(FFLAGS) dsimthnf.f

run: $(OBJS)
	$(FC) $(OBJS)

clean:
	rm -f *.o *.mod run
