##
## Makefile for LIBAOS
##


OFLAGS = -O3
ENERGY_OFLAGS = -O2

LDFLAGS = -lm

SONAME=libaos.so.1.0.2
LIBOBJ = SW_main.o sequence.o SmithWaterman.o input.o energy.o energy_cal.o data.o LZW_Pick70.o
LIBSRC = SW_main.c sequence.c SmithWaterman.c input.c energy.c energy_cal.c data.c LZW_Pick70.c

libaos : $(SONAME) $(LIBOBJ)

$(SONAME) : $(LIBOBJ)
	cc -shared -o $(SONAME) $(LIBOBJ)

$(LIBOBJ) :  $(LIBSRC)
	cc -O2 -fPIC -c $(LIBSRC)

clean :
	/bin/rm -f $(LIBOBJ) $(SONAME)


