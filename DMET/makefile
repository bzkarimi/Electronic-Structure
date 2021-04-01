FC = gfortran


EXE = DMET


SRC = \
	array.f90   \
	allocate-array.f90  \
	integrals.f90   \
	geometry.f90  \
	rhf.f90  \
	fci.f90  \
	bfgs_oep.f90  \
	other-subroutine.f90   \
	DMET.f90


LIBS = -L$(LIB) -llapack -lgfortran -lrefblas

LIB = /Users/Borna/Desktop/Graduate/Jason\ Goodpaster/Lapack/lapack-3.6.0/

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=    $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:	
	rm -f $(OBJ) $(EXE) core 
