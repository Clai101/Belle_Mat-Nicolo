DEBUG=-O3

ROOT_INC=$(shell root-config --cflags)
ROOT_LIB=$(shell root-config --libs)

INC=-I/usr/include $(ROOT_INC) -I. 
LIBS=-lm -lstdc++ $(ROOT_LIB) -L/usr/lib64

OPTFLAGS= -pipe -Wall -fPIC -std=c++11

all: test.exe 

test.exe: test.o 
	gcc $(OPTFLAGS)  $(LIBS)  $^ -o  $@


%.o: %.cc
	gcc -c $(OPTFLAGS) $(DEBUG) -MMD $(INC) -o $@ $<

clean:
	-rm *.o *.d *.so *.a test*.root *.exe *_dict.cc
.SUFFIXES:

include $(wildcard *.d)
