#Makefile Projet HPC
#Compilateur MPI
CC = mpicc

CFLAGS = -Wall -O3 -fopenmp

TARGET = mitm

#Fichier source
SRC = mitm.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)


clean:
	rm -f $(TARGET)
	rm -f *.o

.PHONY: all clean