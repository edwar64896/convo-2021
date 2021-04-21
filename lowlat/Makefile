CC 		= gcc -g
CFLAGS 	= -I. $(DEFS)
#DEFS 	= -D _DBG_CONVO
DEFS 	= 
DEPS 	= pipe.h
OBJ 	= pipe.o main.o direct.o

LIBS=-lsndfile -lfftw3_threads -lfftw3 -lm -lpthread -lportaudio

all: llat

clean:
	rm -f *.o llat

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

llat: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

%.raw: $.wav
	sndfile-convert 

