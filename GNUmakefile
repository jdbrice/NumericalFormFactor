OBJS = NumericalFormFactors.o
EXE = nff

#OBJS = MuDst.o createphiwgt.o
#EXE = createphiwgt

#OBJS = MuDst.o StPhysicalHelixD.o StHelixD.o StTofGeometry.o createphiwgt.o
#EXE = createphiwgt

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++
FLAGS = -std=c++11 -fno-inline -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -fno-inline -Wl -lpthread -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)
%.o: %.cpp
	$(COMPILE)  $< 

