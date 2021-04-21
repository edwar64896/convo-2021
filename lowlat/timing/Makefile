CC 		= gcc -g
CFLAGS 	= -I. $(DEFS)
#DEFS 	= -D _DBG_CONVO
DEFS 	= 
DEPS 	= 
OBJ 	= main.o direct.o

LIBS=-lsndfile -lfftw3_threads -lfftw3 -lm -lpthread -lportaudio

all: tcheck

clean:
	rm -f *.o tcheck

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

tcheck: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

%.raw: $.wav
	sndfile-convert 

