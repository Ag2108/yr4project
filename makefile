# Fortran compiler
FC = gfortran

# Compiler flags
FFLAGS = -O0 -g -Wall -fcheck=all -std=f2008 -fmax-errors=1

# Common module source files
MODULES = Basics.f90 Htn.f90 TightBinding.f90
TEST 	= Test.f90
MODULE_OBJS = $(MODULES:.f90=.o)
INPUT="read.in"

# Main program
MAIN1 = Main.f90
MAIN1_OBJ = $(MAIN1:.f90=.o)

MAIN2 = Test.f90
MAIN2_OBJ = $(MAIN2:.f90=.o)

# Executable name
EXE1 = mytb
EXE2 = exe_test

# Default build target (only build what exists)
all: $(EXE1)

# Build File1
$(EXE1): $(MODULE_OBJS) $(MAIN1_OBJ)
	$(FC) $(FFLAGS) -o $@ $(MODULE_OBJS) $(MAIN1_OBJ) -llapack -lblas

$(EXE2): $(MODULE_OBJS) $(MAIN2_OBJ)
	$(FC) $(FFLAGS) -o $@ $(MODULE_OBJS) $(MAIN2_OBJ) -llapack -lblas

# Compile each module or source file to object file
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

run: $(EXE1)
	echo $(INPUT) | ./$(EXE1) > out.txt

comp: $(EXE1)
test: $(EXE2)
	./$(EXE2)

# Clean build artifacts
clean:
	rm -f *.o *.mod $(EXE1) $(EXE2)

.PHONY: all run comp test clean
