#
CC        = g++
CLINKER   = $(CC)
OPTFLAGS  = 
LIBS      = -lfftw3 -lm
CFLAGS    = $(OPTFLAGS)
LINKFLAGS = $(OPTFLAGS) $(LIBS) 

all: fftest 1dtest

fftest: fft_test.o
	$(CLINKER) -o fftest fft_test.o $(LINKFLAGS) 

1dtest: 1d_fftw_test.o
	$(CLINKER) -o 1dtest 1d_fftw_test.o $(LINKFLAGS) 

clean:
	rm -r fftest 1dtest *.o
.c.o:
	$(CC) $(CFLAGS) -c $<

.SUFFIXES:	.c .C
