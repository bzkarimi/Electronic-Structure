FC = gfortran


EXE = bfgs


SRC = bfgs.f

LIBS = -L$(LIB) -llapack -lgfortran -lrefblas

LIB = /Users/Borna/Desktop/Graduate/Jason\ Goodpaster/Lapack/lapack-3.6.0/

.SUFFIXES:
.SUFFIXES: .f .o

OBJ=    $(SRC:.f=.o)

.f.o:
	$(FC) $(FFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:	
	rm -f $(OBJ) $(EXE) core 
