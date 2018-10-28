CC=g++
CFLAGS=-c -Wall -std=c++11 -O3 -g
INCLUDES=-I. -I.. -I$(shell root-config --incdir)
#(root-config --incdir)
LIBDIRS=-L/usr/lib
LIBS=-L../lib -lRhoNNO
LDFLAGS=$(shell root-config --glibs)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

SOURCES1=ReadEvent.cxx Geo.cxx DataStructures.cxx 
SOURCES=main.cxx
SHLIB=Tracker.cxx Point.cxx XMLP.cxx

HEADERS=$(SRC:.cxx=.h)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=tracker
LIBRARY=libprocess.so

all: $(SOURCES) $(LIBRARY) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(LIBRARY)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBRARY) $(ROOTLIBS) $(LIBS) -o $@

$(LIBRARY): Tracker.o XMLP.o XMLP.h Point.h Tracker.h
	clang++ -std=c++11 -fPIC -shared -o $(LIBRARY) XMLP.cxx Point.cxx Tracker.cxx

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm reco $(OBJECTS)

cleanall:
	rm reco $(OBJECTS) *.~ *.root
