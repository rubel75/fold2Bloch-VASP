# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags. Keep '-assume byterecl' !!!
FC = ifort
FLFLAGS = # none
FCFLAGS = -free -assume byterecl # performance
#FCFLAGS = -free -assume byterecl -g -traceback -check all -debug all # debug

# gfortran compiler
#FC = gfortran
#FLFLAGS = # none
#FCFLAGS = -ffree-form # performance
#FCFLAGS = -ffree-form -g -fbacktrace -Wconversion -fcheck=all # debug

# ~~~ Do not edit after that line ~~~

PROGRAM = fold2Bloch

# source files and objects
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90)) \
    $(patsubst %.h, %.mod, $(wildcard *.h))


all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

#%.mod: %.h
#	$(FC) $(FCFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
