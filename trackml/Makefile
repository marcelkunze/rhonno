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
OBJECTS1=ReadEvent.o Geo.o DataStructures.o

SOURCES=main.cxx train.cxx
SHLIB=Tracker.cxx data.cxx topquark.cxx PolarModule.cxx Point.cxx find.cxx swimmer.cxx XMLP.cxx

HEADERS=$(SRC:.cxx=.h)
OBJECTS=$(SOURCES:.cxx=.o) $(SHLIB:.cxx=.o)
EXECUTABLE=tracker
LIBRARY=libtrackml.so

all: $(SOURCES) $(LIBRARY) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(LIBRARY)
	$(CC) -Wl,-stack_size -Wl,0x1000000 $(LDFLAGS) $(OBJECTS) $(LIBRARY) $(ROOTLIBS) $(LIBS) -o $@

ReadEvent: $(OBJECTS1) $(LIBRARY)
	$(CC) $(LDFLAGS) $(OBJECTS1) $(LIBRARY) $(ROOTLIBS) $(LIBS) -o $@

$(LIBRARY): Tracker.o data.o topquark.o PolarModule.o Point.o find.o swimmer.o XMLP.o XMLP.h Point.h Tracker.h
	clang++ -std=c++11 -fPIC -shared -o $(LIBRARY) $(SHLIB)

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm Readevent $(EXECUTABLE) $(LIBRARY) $(OBJECTS) $(OBJECTS1)
