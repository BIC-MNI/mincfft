PROGS = mincfft
OBJS = $(PROGS:=.o) fft_support.o

CC=cc

OPTIONS = -O3
INCLUDES = -I/usr/local/include
CFLAGS = $(OPTIONS) $(INCLUDES)

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lfftw -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS)


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
