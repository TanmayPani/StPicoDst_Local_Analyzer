# Define compiler
CXX := g++

ROOTCFLAGS	    := $(shell root-config --cflags)
PICOPATH        := /Users/tanmaypani/StPicoEvent
PICOCFLAGS      := -I$(PICOPATH)
FJCFLAGS        := $(shell fastjet-config --cxxflags)
FJWRAPPERPATH   := /Users/tanmaypani/FASTJET/FastJetWrapper
FJWRAPPERCFLAGS := -I$(FJWRAPPERPATH)

ROOTLIBS        := $(shell root-config --libs)
PICOLIBS        = -Wl,-rpath,$(PICOPATH) -L$(PICOPATH) -lStPicoDst
FJLIBS	        := $(shell fastjet-config --libs)
FJWRAPPERLIBS   = -Wl,-rpath,$(FJWRAPPERPATH) -L$(FJWRAPPERPATH) -lFastJetWrapper

ROOTINC         := $(shell root-config --incdir)

CFLAGS = $(ROOTCFLAGS) -I. $(PICOCFLAGS) $(FJCFLAGS) $(FJWRAPPERCFLAGS) -O2 -fPIC -Wall -W -Woverloaded-virtual -Wno-deprecated-declarations
CFLAGS += -pipe -std=c++14 -D_VANILLA_ROOT_ 

LIBS = $(ROOTLIBS) $(PICOLIBS) $(FJLIBS) $(FJWRAPPERLIBS)

INCS = $(ROOTINC) $(PICOCFLAGS) $(FJCFLAGS) $(FJWRAPPERCFLAGS)

# Define output library
STPICOANALIB := libStPicoAnalyzer.dylib

# Compile all *.cpp classes in the directory
SRC := $(shell find . -name "*.cpp")

all: $(STPICOANALIB)

# $(SRC:.cc=.o)
$(STPICOANALIB): $(SRC:.cpp=.o) StPicoAnalyzer_Dict.C
	$(CXX) $(CFLAGS) -shared $^ -o $(STPICOANALIB) $(LIBS) 

%.o: %.cpp
	$(CXX) -fPIC $(CFLAGS) -c -o $@ $<

StPicoAnalyzer_Dict.C: $(shell find . -name "*.h" ! -name "*LinkDef*")
	rootcint -f $@ -c -D_VANILLA_ROOT_ -DROOT_CINT -D__ROOT__ -I. -I$(INCS) $^ StPicoAnalyzer_LinkDef.h

.PHONY: clean distclean

clean:
	rm -vf *.o StPicoAnalyzer_Dict*

distclean:
	rm -vf *.o StPicoAnalyzer_Dict* $(STPICOANALIB)
check:
	@echo "CXX = $(CXX)"
	@echo "CFLAGS = $(CFLAGS)"
	@echo "LIBS = $(LIBS)"
	@echo "INCS = $(INCS)"
	@echo "SRC = $(SRC)"
	
