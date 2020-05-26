# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags. Keep '-assume byterecl' !!!
FC = ifort
FCFLAGS = -assume byterecl #-g -traceback -check all -debug all
FLFLAGS = -assume byterecl #-g -traceback -check all -debug all

# ~~~ Do not edit after that line ~~~

PROGRAM = fold2Bloch

# source files and objects
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90)) \
    $(patsubst %.h, %.mod, $(wildcard *.h))


all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FLFLAGS) -c $<

#%.mod: %.h
#	$(FC) $(FLFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
