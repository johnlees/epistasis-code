# Default prefix. Can change this to system directories if needed
PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

# Intel compiler - uncomment if you have icpc and mkl.
CXX=icpc
CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
SEER_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -mkl
# gcc
#CXXFLAGS=-Wall -O3 -std=c++11
#SEER_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -llapack -lblas -lpthread

CPPFLAGS=-I$(PREFIX)/include -I../gzstream -I../dlib -I/usr/local/hdf5/include -D DLIB_NO_GUI_SUPPORT=1 -D DLIB_USE_BLAS=1 -D DLIB_USE_LAPACK=1 -DARMA_USE_HDF5=1

PROGRAMS=epistasis

OBJECTS=pair.o stats.o logitFunction.o logisticRegression.o common.o cmdLine.o epistasis.o

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

epistasis: $(OBJECTS)
	$(LINK.cpp) $^ $(SEER_LDLIBS) -o $@

.PHONY: all clean install

