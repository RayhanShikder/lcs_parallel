# Sample OpenMP Makefile

###############################################
# Change these files names to match your own: #
###############################################

# name of the executable to be generated
PROG = row_wise_v1
# space-delimited list of header files
HDRS =
# space-delimited list of source files
SRCS = row_wise_v1.c 

#######################
# Don't change these: #
#######################

# directory to store object files
OBJDIR = object
# names of object files
OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(SRCS)) 

# name of the compiler
CC = gcc
# additional compiler flags to pass in
CFLAGS = -fopenmp -Wall --std=c99 -L.
# libraries for the linker
LIBS = -lm

####################
# Compiling rules: #
####################
# WARNING: *must* have a tab before each definition

# invoked when "make" is run
all : $(OBJDIR) $(PROG)

# links object files into executable
$(PROG) : $(OBJS)
	$(CC) $(CFLAGS) $^ -o $(PROG) $(LIBS)

# compiles source files into object files
object/%.o : %.c $(HDRS)
	$(CC) -c $(CFLAGS) $< -o $@ $(LIBS)

# creates directory to store object files
$(OBJDIR) :
	mkdir -p $@/

# cleans up object files and executable
# type "make clean" to use
# Note: you can add your own commands to remove other things (e.g. output files)
clean:
	rm -rf object/
	rm -rf $(PROG)
	rm -rf output.txt
	rm -rf Q1.o
