OBJS = StPhysicalHelix.o StHelix.o Lcp_Reconstruction.o
EXE = Lcp_Reconstruction

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(shell root-config --incdir)
LDFLAGS = -L$(shell root-config --libdir)

CXX = g++
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $<

clean:
	rm $(OBJS)

