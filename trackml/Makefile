CC=g++
CFLAGS=-c -Wall -std=c++11 -O3 -g
INCLUDES=-I. -I.. -I$(shell root-config --incdir)
#(root-config --incdir)
LIBDIRS=-L/usr/lib
LIBS=-L../lib -lRhoNNO
LDFLAGS=$(shell root-config --glibs)

#SOURCES=ReadEvent.cxx Geo.cxx DataStructures.cxx xmlp.cxx netbase.cxx percept.cxx net.cxx
SOURCES=ReadEvent.cxx Geo.cxx DataStructures.cxx 
SHLIB=findTracks.cxx XMLP.cxx

HEADERS=$(SRC:.cxx=.h)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=ReadEvent
LIBRARY=libprocess.so

all: $(SOURCES) $(LIBRARY) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBRARY) $(LIBS) -o $@

$(LIBRARY): findtracks.o XMLP.o
	clang++ -fPIC -shared -o $(LIBRARY) XMLP.cxx findtracks.cxx

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm reco $(OBJECTS)

cleanall:
	rm reco $(OBJECTS) *.~ *.root
