# Makefile for g3MPI test

CSRC= g3MPI.c main.c g3transpose.c
OBJ=$(CSRC:.c=.o)

HSRC= grid.h mytypes.h

.SUFFIXES:
.SUFFIXES: .c .o


.c.o:
	$(CC) $(CFLAGS) -c $<

fftpf3d: $(OBJ)
	$(CC) $(LDFLAGS) -o fftpf3d -DORIG_ORDER $(OBJ) -lm $(LIBS)
	-rm -f $(OBJ)

clean:
	-rm -f a.out $(OBJ) *core

$(OBJ): $(HSRC)
