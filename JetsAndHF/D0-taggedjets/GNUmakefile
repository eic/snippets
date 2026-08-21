OBJS = StPhysicalHelix.o StHelix.o analysis.o
EXE = analysis

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(shell root-config --incdir)
LDFLAGS = -L$(shell root-config --libdir)

FASTJETCFLAGS   = $(shell fastjet-config --cxxflags)
FASTJETLIBS     = $(shell fastjet-config --libs)


CXX = g++-14
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS) $(FASTJETCFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS) $(FASTJETLIBS)

%.o: %.cxx
	$(COMPILE)  $<

clean:
	rm $(OBJS)

