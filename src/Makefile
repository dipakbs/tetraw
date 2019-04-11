CC=gcc
PROFFLAGS=-pg
DEBUG=-O0 -g -DDEBUG 
MAIN=tetrahedral
LDFLAGS=-lm

tetrahedral : tetrahedral.c
	$(CC) $(LDFLAGS) -O2 -o $(MAIN) tetrahedral.c 

profile:  tetrahedral.c
	$(CC) $(PROFFLAGS)  $(LDFLAGS) tetrahedral.c
debug:	tetrahedral.c
	$(CC) $(DEBUG) $(LDFLAGS) -o $(MAIN) tetrahedral.c
clean: 
	rm $(MAIN)
