# Compiler
FC = gfortran

# Compiler flags
# Added -cpp to run the C preprocessor which is required for #include statements in .f90 files
FFLAGS = -O2 -g -Wall -I.. 

# Paths to SuperLU
SUPERLU_DIR = ../superlu/build
SUPERLU_INC = -I$(SUPERLU_DIR)/SRC
SUPERLULIB   	= $(SUPERLU_DIR)/SRC/libsuperlu.a

LIBS = $(SUPERLULIB) /opt/homebrew/opt/openblas/lib/libopenblas.dylib -lm
# Source files
SRCS = constants.f90 interpolation.f90 gas.f90 dust.f90 tripod.f90
MAIN = main_test.f90

#SuperLU object files (you may need to adjust this based on your SuperLU build)
splu_objs = $(SUPERLU_DIR)/../FORTRAN/c_fortran_dgssv.o

# Object files automatically generated from SRCS
OBJS = $(SRCS:.f90=.o)
MAINOBJ = $(MAIN:.f90=.o)

# Module files generated
MODS = $(SRCS:.f90=.mod)

# Target library (or change to an executable if you add a main program)
TARGET = libtripod.a

# Default target
all: $(TARGET)

# Create a static library
$(TARGET): $(OBJS) $(splu_objs)
	ar rcs $@ $(OBJS) $(splu_objs)

# Compile Fortran files
%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

# Dependencies (if modules depend on each other, define them here)
# Example: tripod.o: dust.o gas.o

# All object files depend on parameters.h being up to date
$(OBJS): ../parameters.h

dust.o: constants.o
gas.o: interpolation.o
tripod.o: dust.o gas.o ${splu_objs}

main: $(MAINOBJ) $(TARGET)
	$(FC) $(FFLAGS) -o $@ $(MAINOBJ) $(OBJS) libtripod.a  $(splu_objs) $(LIBS)

clean:
	rm -f $(OBJS) $(MODS) $(TARGET)

.PHONY: all clean
