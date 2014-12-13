
CC = gcc
CFLAGS = -O3 -g 
LIBS = -lm -lgsl -lgslcblas


HEADERS = mg_func.h
OBJECTS = mg_func.o main.o 
PROGRAM = mg

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c


