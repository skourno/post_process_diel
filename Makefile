
.SUFFIXES:	(.SUFFIXES) .f90 .f

## Compilers and flags ##
FC      = gfortran
FCFLAGS =  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

## Linker and flags ##
LD      = $(FC)
LDFLAGS =  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
LIBS    =

## Compilation directives ##

.f90.o:
	$(FC) -c $(FCFLAGS) $*.f90

MOBJECTS = main.o

BOBJECTS=nrtype.o\
	system.o\
	utilities.o\
	sph_harmonics.o\
	read_cases.o\
	input.o\
	histograms.o\
	g_r.o\
	compute.o\
	flow.o

OBJECTS =  $(BOBJECTS) $(MOBJECTS)  

## Executable ##
CMD =  diel

## Make directives ##

$(CMD):	$(OBJECTS)
	$(LD) -o $(CMD) $(LDFLAGS) $(OBJECTS) $(LIBS)

clean:
	/bin/rm *.o *.mod

distclean:
	/bin/rm *.o *.mod $(CMD)

