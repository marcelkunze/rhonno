CC=g++
CFLAGS=-c -Wall -std=c++11 -O3 -g
INCLUDES=-I. -I.. -I$(shell root-config --incdir)
#(root-config --incdir)
LIBDIRS=-L/usr/lib
LIBS=-L../lib -lRhoNNO
LDFLAGS=$(shell root-config --glibs)

#SOURCES=ReadEvent.cxx Geo.cxx DataStructures.cxx xmlp.cxx netbase.cxx percept.cxx net.cxx
SOURCES=ReadEvent.cxx Geo.cxx DataStructures.cxx

HEADERS=$(SRC:.cxx=.h)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=ReadEvent

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm reco $(OBJECTS)

cleanall:
	rm reco $(OBJECTS) *.~ *.root
