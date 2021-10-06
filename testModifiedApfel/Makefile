
CXX = g++

CXXFLAGS += -O3 -fPIC -std=c++11

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELPPINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELPPINCS)

CLIBS += $(LHAPDFLIBS) $(APFELPPLIBS)

install : all
all : Evolution StructureFunctions

Evolution: Evolution.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

StructureFunctions: StructureFunctions.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la Evolution StructureFunctions *~

