CC 		= gcc -g
CFLAGS 	= -I. $(DEFS)
#DEFS 	= -D _DBG_CONVO
DEFS 	= 
DEPS 	= convo.h
OBJ 	= convo.o  oacon.o 

LIBS=-lsndfile -lfftw3_threads -lfftw3 -lm

all: convo

clean:
	rm -f *.o dbg* output* convo iconvo impulseresponse*

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

convo: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

%.raw: $.wav
	sndfile-convert 

