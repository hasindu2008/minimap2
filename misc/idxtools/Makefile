CC       = gcc
CFLAGS   = -g -Wall -O2

SRC = divide.c
BINARY = divide

.PHONY: clean test

$(BINARY) : $(SRC) 
	$(CC) $(CFLAGS) $(SRC) -o $@ 
	
clean: 
	rm -rf divide *.o *.out

test: $(BINARY)
	cd tmp && ./test.sh 
