# Makefile for Linux
# Assumes the environment variable ROOTSYS to exist
#

# Define directory paths

ifeq ($(RHO),)
RHO = ..
endif

OBJDIR        = $(RHO)/tmp
LIBDIR        = $(RHO)/lib
BINDIR        = $(RHO)/bin

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
# Linux with gcc
CXX           = c++
CXXFLAGS      = -I$(RHO) -O2 -Wall -fPIC -pthread -std=c++11 -m64
LD            = c++
LDFLAGS       = -Wl,-rpath,/lib:$(LIBDIR):$(shell root-config --libdir)
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -L$(LIBDIR) -lRhoNNO -lTreePlayer
GLIBS         = $(ROOTGLIBS)

# Get dependencies:
include GNUmake-depend.mk

strain:
	$(CXX) $(CXXFLAGS) strain.cxx -o $(RHO)/bin/strain $(LDFLAGS) $(LIBS)

utrain:
	$(CXX) $(CXXFLAGS) utrain.cxx -o $(RHO)/bin/utrain $(LDFLAGS) $(LIBS)

NetworkTrainer:
	$(CXX) $(CXXFLAGS) NetworkTrainer.cxx -o $(RHO)/bin/NetworkTrainer $(LDFLAGS) $(LIBS)

