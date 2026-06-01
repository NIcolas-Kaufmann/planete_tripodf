# Compiler
FC = gfortran

# Compiler flags
# Added -cpp to run the C preprocessor which is required for #include statements in .f90 files
FFLAGS = -O2 -g -Wall -I.. 

# Paths to SuperLU
SUPERLU_DIR = ../superlu
SUPERLU_INC = -I$(SUPERLU_DIR)/SRC
SUPERLU_LIB = -L$(SUPERLU_DIR)/build/SRC -lsuperlu

# LAPACK/BLAS are required by SuperLU
BLAS_LIB = -lblas -llapack

# Source files
SRCS = constants.f90 interpolation.f90 gas.f90 dust.f90 tripod.f90

# Object files automatically generated from SRCS
OBJS = $(SRCS:.f90=.o)

# Module files generated
MODS = $(SRCS:.f90=.mod)

# Target library (or change to an executable if you add a main program)
TARGET = libtripod.a

# Default target
all: $(TARGET)

# Create a static library
$(TARGET): $(OBJS)
	ar rcs $@ $(OBJS)

# Compile Fortran files
%.o: %.f90
	$(FC) $(FFLAGS) $(SUPERLU_INC) -c $< -o $@

# Dependencies (if modules depend on each other, define them here)
# Example: tripod.o: dust.o gas.o

# All object files depend on parameters.h being up to date
$(OBJS): ../parameters.h

dust.o: constants.o
gas.o: interpolation.o
tripod.o: dust.o gas.o

clean:
	rm -f $(OBJS) $(MODS) $(TARGET)

.PHONY: all clean
