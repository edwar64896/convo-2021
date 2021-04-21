CC 		= gcc -g
CFLAGS 	= -I. $(DEFS)
DEFS 	= -D _DBG_CONVO
DEPS 	= convo.h
OBJ 	= convo.o  oacon.o 
OBJ2 	= ioacon.o iconvo.o

LIBS=-lsndfile -lfftw3_threads -lfftw3 -lm

all: convo iconvo

clean:
	rm -f *.o dbg* output* convo iconvo impulseresponse*

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

convo: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

iconvo: $(OBJ2)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

%.raw: $.wav
	sndfile-convert 

