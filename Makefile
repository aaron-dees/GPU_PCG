CC = gcc
CFLAGS=-std=gnu99 -O0 -march=native -Wall -m64 -g


pcgSolver: main.o params.o matvec_mul.o pcg_algo.o precondit.o
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS) -pg -lm

.c.o:
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -pg -lm

.PHONY: clean
clean:
	rm -f *.o *~ pcgSolver
