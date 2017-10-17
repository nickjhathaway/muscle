# Porting notes:
# For non-Intel platforms, change the -march option as appropriate.
# For Solaris and other platforms where the logf function
# is missing from the math library, add the following line
# to the end of muscle.h:
# #define logf(x)	((float) log(x))
# Using -static increases the executable size and thus gives a very
# small increase in start time, but is more portable (the binding
# to dynamic libraries often breaks when a new library is released).

CFLAGS = -O3 -march=pentiumpro -mcpu=pentiumpro -funroll-loops -Winline -DNDEBUG=1
LDLIBS = -lm -static
# LDLIBS = -lm

OBJ = .o
EXE =

RM = rm -f
CP = cp

GPP = g++
LD = $(GPP) $(CFLAGS)
CPP = $(GPP) -c $(CFLAGS) 
CC = gcc -c $(CFLAGS) 

all: muscle

CPPSRC = $(sort $(wildcard *.cpp))
CPPOBJ	= $(subst .cpp,.o,$(CPPSRC))

$(CPPOBJ): %.o: %.cpp
	$(CPP) $< -o $@

muscle: $(CPPOBJ)
	$(LD) -o muscle $(CPPOBJ) $(LDLIBS)
	strip muscle
