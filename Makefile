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
