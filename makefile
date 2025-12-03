# Fortran compiler
FC = gfortran

# Compiler flags
FFLAGS = -O0 -g -Wall -fcheck=all -std=f2008 -fmax-errors=1

# Common module source files
MODULES = Basics.f90 Htn.f90 TightBinding.f90
MODULE_OBJS = $(MODULES:.f90=.o)
INPUT="read.in"

# Main program
MAIN1 = Main.f90
MAIN1_OBJ = $(MAIN1:.f90=.o)

# Executable name
EXE1 = mytb

# Default build target (only build what exists)
all: $(EXE1)

# Build File1
$(EXE1): $(MODULE_OBJS) $(MAIN1_OBJ)
	$(FC) $(FFLAGS) -o $@ $(MODULE_OBJS) $(MAIN1_OBJ) -llapack -lblas

# Compile each module or source file to object file
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Optional: pass command-line args to Task1: make run ARGS="a b c d"
ARGS ?=
run: $(EXE1)
	echo $(INPUT) | ./$(EXE1) $(ARGS)

# Keep run1 as an alias if you like
run1: run

# Clean build artifacts
clean:
	rm -f *.o *.mod $(EXE1)

.PHONY: all run run1 clean
