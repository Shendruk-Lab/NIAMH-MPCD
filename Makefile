#===================================================================================================
# Make file for the MPC simulation program
# Tyler Shendruk
#====================================================================================================

defines		= -DPHASE_COUNT=2 -DTYPE_COUNT=3 -DPOLY_SETS=1 -DQ_SETS=4
conditional = -DTHERMOSTAT_DPD

# name of program
program    = mpcd/mpcd.out
progName   = $(program)

# shell commands to automatically build list of files from *.c and *.h in directory
sources    = $(shell ls mpcd/subroutines/*.c) $(program:.out=.c) $(shell ls md/*.c) $(shell ls dependencies/cJson/*.c)
headers    = $(shell ls mpcd/headers/*.h) $(shell ls md/*.h) $(shell ls dependencies/cJson/*.h)
objects    = $(sources:.c=.o)

# other variables
cc         :=  gcc
cflags     :=  -Wall
lib        :=  -lm
opt        :=  -O3

#----------------------------------------------------------------------------------------------------
# By default, make will build the first target, $(program) here.
# If any of the $(objects) are more recent than $(program), then linking takes place
$(program): $(objects)
	@echo "LINKING:  $(program)"
	@$(cc) $(cflags) $(opt) $(objects) $(defines) $(conditional) $(lib) -o $(progName)
	mv $(progName) ./
#----------------------------------------------------------------------------------------------------
# If any .c file is more recent than its .o file, then they are compiled
%.o:    %.c
	@echo "COMPILING: $(cc) $(cflags) $(opt) $<"
	@$(cc) -c $(cflags) $(defines) $(conditional) $(opt)  $< -o $@
#----------------------------------------------------------------------------------------------------
# Use the .PHONY declaration when you want to prevent the target ("debug" here) from being
# interpreted as the name of a file.
.PHONY:    debug
debug:
	make -e opt="-g" progName="mpcd/mpcdDebug.out"
#----------------------------------------------------------------------------------------------------
# phony for extra debugging/ valgrind
.PHONY:    debug+
debug+:
	make -e opt="-ggdb3" progName="mpcd/mpcdDebugPlus.out"
#----------------------------------------------------------------------------------------------------
# phony for profiling
.PHONY:    prof
prof:
	make -e opt="$(opt) -pg" progName="mpcd/mpcdProf.out"
#----------------------------------------------------------------------------------------------------
# phony for OMP
.PHONY:    omp
prof:
	make -e lib="$(lib) -fopenmp" progName="mpcd/mpcdOMP.out"
#----------------------------------------------------------------------------------------------------
.PHONY: clean
clean:
	@echo "Removing executable and object files"
# 	@echo "	About to perform:  /bin/rm -f $(program) $(objects) "
	@( /bin/rm -f $(program) $(objects) ./*.out )
#----------------------------------------------------------------------------------------------------
.PHONY: docs
docs:    $(sources)
	@echo "Building documentation"
	doxygen doxyfile
#----------------------------------------------------------------------------------------------------
