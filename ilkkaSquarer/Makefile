
BIN     =  ilkkaSquarer

OBJ     =  math.o \
           matrixsquaring.o \
           periodicbc.o \
           init.o \
	   main.o

LIBS = 

FFLAGS = -O2
#FFLAGS = -CB -warn all
FC = gfortran



$(BIN): $(OBJ)
	$(FC) $(FFLAGS) -o $(BIN) $(OBJ) $(LIBS)


%.o : %.f90
	$(FC) -c $(FFLAGS) $<

#.o: 
#	$(FC) $(FLAGS) -c $(.SOURCE)

clean: 
	rm -f core
	rm -f *.mod
	rm -f $(OBJ)
	rm -f $(BIN)
