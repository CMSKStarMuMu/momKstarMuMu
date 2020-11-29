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
ROOFITFLAGS1  := $(shell root-config --cflags --libs) -lMinuit -lRooFit -lRooFitCore -lFoam

Moment : $(EXECUTABLE).cc
	$(CXX) $(CXXFLAG) -o $@  $(ROOTFLAGS) $(ROOFITFLAGS1) -D"ROOFIT" $@.cc

LOCALCC := $(wildcard *.cc)
clean:
	$(RM) $(LOCALCC:%.cc=%)  $(CLASSDICT2)  Macros_C.*
