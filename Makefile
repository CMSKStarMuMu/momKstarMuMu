OS            := $(shell uname)

#ifeq ($(OS), Linux)
CXX           := g++
DEBUG         := -g
#endif

CXXFLAGS      := $(DEBUG) -ansi -Wall -Wextra -m64 -O3 -std=c++11
#DEBUGFLAGS := -O3 -Wall -std=c++0x
#CXXFLAGS := $(DEBUGFLAGS)

EXECUTABLE    := Moment

ROOTFLAGS     := $(shell root-config --glibs)
ROOFITFLAGS   := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore

Moment : $(EXECUTABLE).cc
	$(CXX) $(CXXFLAG) -o $@  $(ROOTFLAGS) $(ROOFITFLAGS1) -D"ROOFIT" $@.cc

LOCALCC := $(wildcard *.cc)
clean:
	$(RM) $(LOCALCC:%.cc=%)  $(CLASSDICT2)  Macros_C.*
