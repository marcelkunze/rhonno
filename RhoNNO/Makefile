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
CXXFLAGS      = $(RHOARCH) -DANSISTREAM -fno-exceptions -fPIC -I.. -pthread -stdlib=libc++ -std=c++11 -m64
LD            = g++
LDFLAGS       = $(RHOARCH) -Wl,-rpath,/lib:$(LIBDIR):$(ROOTSYS)/lib
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lEG -lProof -lTreePlayer -lThread
GLIBS         = $(ROOTGLIBS) -lEG -lProof -lTreePlayer -lThread

# Get dependencies:
include GNUmake-depend.mk

strain:
	$(CXX) -pthread -stdlib=libc++ -std=c++11 -m64 -o $(RHO)/bin/strain -I$(RHO) -I$(ROOTSYS)/include $(LIBS) $(RHO)/lib/libRhoNNO.so strain.cxx

utrain:
	$(CXX) -pthread -stdlib=libc++ -std=c++11 -m64 -o $(RHO)/bin/utrain -I$(RHO) -I$(ROOTSYS)/include $(LIBS) $(RHO)/lib/libRhoNNO.so utrain.cxx

NetworkTrainer:
	$(CXX) -pthread -stdlib=libc++ -std=c++11 -m64 -o $(RHO)/bin/NetworkTrainer -I$(RHO) -I$(ROOTSYS)/include $(LIBS) $(RHO)/lib/libRhoNNO.so NetworkTrainer.cxx

