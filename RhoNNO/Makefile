# Makefile for Linux
# Assumes the environment variable ROOTSYS to exist
#

# Define directory paths

OBJDIR        = $(RHO)/tmp
LIBDIR        = $(RHO)/lib
BINDIR        = $(RHO)/bin

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
# Linux with egcs
CXX           = g++
CXXFLAGS      = $(RHOARCH) -DANSISTREAM -fno-exceptions -fPIC -I..
LD            = g++
LDFLAGS       = $(RHOARCH) -Wl,-rpath,/lib:$(LIBDIR):$(ROOTSYS)/lib
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lEG -lProof -lTreePlayer -lThread
GLIBS         = $(ROOTGLIBS) -lEG -lProof -lTreePlayer -lThread

# Get dependencies:
include GNUmake-depend.mk

strain:
	$(CXX) -o $(RHO)/bin/strain -I$(RHO) -I$(ROOTSYS)/include $(LIBS) $(RHO)/lib/libRhoNNO.so strain.cxx

NetworkTrainer:
	$(CXX) -o $(RHO)/bin/NetworkTrainer -I$(RHO) -I$(ROOTSYS)/include $(LIBS) $(RHO)/lib/libRhoNNO.so NetworkTrainer.cxx

