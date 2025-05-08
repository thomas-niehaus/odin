# Fortran compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall

# Source directory
SRC_DIR = src

# Object files
OBJS = $(SRC_DIR)/odin.o \
       $(SRC_DIR)/gettab.o \
       $(SRC_DIR)/slkode.o \
       $(SRC_DIR)/skspar.o \
       $(SRC_DIR)/skss.o \
       $(SRC_DIR)/sksp.o \
       $(SRC_DIR)/sksd.o \
       $(SRC_DIR)/skpp.o \
       $(SRC_DIR)/skpd.o \
       $(SRC_DIR)/skdd.o \
       $(SRC_DIR)/selfs.o \
       $(SRC_DIR)/selfp.o \
       $(SRC_DIR)/selfd.o

# Executable name
EXE = odin

# Installation prefix (can be overridden when calling make)
PREFIX ?= /usr/local/bin

# Default target
all: $(EXE)

# Build executable
$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Compile source files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod $(EXE)

# Install rule
install: $(EXE)
	install -d $(PREFIX)
	install -m 755 $(EXE) $(PREFIX)

.PHONY: all clean install
