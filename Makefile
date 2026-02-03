# Makefile for Linux
# Assumes the environment variable ROOTSYS to exist
#

# Define directory paths

ifeq ($(RHO),)
RHO = .
endif

OBJDIR        = $(RHO)/tmp
LIBDIR        = $(RHO)/lib
BINDIR        = $(RHO)/bin
DATADIR       = $(RHO)/Data

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
# Linux with clang
CXX           = clang++
CXXFLAGS      = -I$(RHO) -O2 -Wall -fPIC -pthread -std=c++17 -m64 -Wunused-parameter -Wshadow
LD            = clang++
LDFLAGS       = -Wl,-rpath,/lib:$(LIBDIR):$(shell root-config --libdir)
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -L$(LIBDIR) -lRhoNNO -lTreePlayer
GLIBS         = $(ROOTGLIBS)

# Get dependencies:
include GNUmake-depend.mk

all: strain utrain NetworkTrainer NNOTracker RadonTracker

testall: teststrain testutrain testNetworkTrainer testNNOTracker testRadonTracker

strain: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) strain.cxx -o $(RHO)/bin/strain $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

teststrain: strain
	. $(PWD)/environment; $(BINDIR)/strain

utrain: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) utrain.cxx -o $(RHO)/bin/utrain $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

testutrain: utrain
	. $(PWD)/environment; $(BINDIR)/utrain $(DATADIR)/ppe.dat

NetworkTrainer: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) NetworkTrainer.cxx -o $(RHO)/bin/NetworkTrainer $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

testNetworkTrainer: NetworkTrainer
	. $(PWD)/environment; $(BINDIR)/NetworkTrainer pid.nno

NNOTracker: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) NNOTracker.cxx -o $(RHO)/bin/NNOTracker $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

testNNOTracker: NNOTracker
	. $(PWD)/environment; $(BINDIR)/NNOTracker $(DATADIR)/event2

RadonTracker: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) RadonTracker.cxx -o $(RHO)/bin/RadonTracker $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

testRadonTracker: RadonTracker
	. $(PWD)/environment; $(BINDIR)/RadonTracker

RadonOptimizer: $(LIBA) $(LIBSO)
	$(CXX) $(CXXFLAGS) RadonOptimizer.cxx -o $(RHO)/bin/RadonOptimizer $(LDFLAGS) $(LIBS)
	mv *.pcm $(LIBDIR)

testRadonOptimizer: RadonOptimizer
	. $(PWD)/environment; $(BINDIR)/RadonOptimizer
