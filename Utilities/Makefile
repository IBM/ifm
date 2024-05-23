## Makefile

.PHONY: clean realclean

### these flags can be modified
#OPTFLAG = -g 
ifeq ("$(shell uname)","Darwin")
OPTFLAG = -g -O3 -fno-inline -I/opt/X11/include -Wno-format
else
OPTFLAG = -g -O3 -fno-inline
endif

### if everything goes fine, should not touch anything below

ifeq ("$(shell uname -m)","i686")
    BITS = -m32
else
    BITS = -m64
endif

CXXFLAGS      = -I. -I.. -I../Core -I../FileIO -Wall ${OPTFLAG} ${OPTIONS} ${BITS}
HDRS          = $(wildcard *.h)
OBJS          = $(patsubst %.C,%.o, $(wildcard *.C))
NETCDF_OBJS   = ../FileIO/netcdf_io.o
NETCDF_LIBS   = -lnetcdf_c++ -lnetcdf
ifeq ("$(shell uname)","Darwin")
LIB_DIRS      = -L/opt/X11/lib
else
LIB_DIRS      = 
endif

PAINTER_OBJS  = painter.o $(NETCDF_OBJS)
PAINTER_LIBS  = -lm -lgd -lpthread $(NETCDF_LIBS) $(LIB_DIRS)

GETVALUE_OBJS = getvalue.o ../FileIO/interpolation.o $(NETCDF_OBJS)
GETVALUE_LIBS = -lm $(NETCDF_LIBS) $(LIB_DIRS)

HISTOGRAM_MAIN_OBJ = histogramMain.o
HISTOGRAM_OBJS = histogram.o ../FileIO/interpolation.o $(NETCDF_OBJS)
HISTOGRAM_LIBS = -lm $(NETCDF_LIBS) $(LIB_DIRS)

DEMVIEWER_OBJS = demviewer.o $(NETCDF_OBJS)
ifeq ("$(shell uname)","Darwin")
DEMVIEWER_LIBS = -lm -lGL -lGLU -lglut -lSDL -lpthread $(NETCDF_LIBS) $(LIB_DIRS) -lSDLmain -lSDL -framework Cocoa
else
DEMVIEWER_LIBS = -lm -lGL -lGLU -lglut -lSDL -lpthread $(NETCDF_LIBS) $(LIB_DIRS)
endif

RISKESTIMATOR_OBJS = riskestimator.o $(NETCDF_OBJS)
RISKESTIMATOR_LIBS = -lm $(NETCDF_LIBS) $(LIB_DIRS)

## change this line since we have different top level drivers
TARGETS = painter getvalue histogram demviewer riskestimator
#TARGETS = painter getvalue histogram riskestimator

all: $(TARGETS)

painter: $(PAINTER_OBJS)
	$(CXX) $^ -o $@ ${OPTFLAG} ${PAINTER_LIBS}

getvalue: $(GETVALUE_OBJS)
	$(CXX) $^ -o $@ ${OPTFLAG} ${GETVALUE_LIBS}

histogram: $(HISTOGRAM_OBJS) $(HISTOGRAM_MAIN_OBJ)
	$(CXX) $^ -o $@ ${OPTFLAG} ${HISTOGRAM_LIBS}

demviewer: $(DEMVIEWER_OBJS)
	$(CXX) $^ -o $@ ${OPTFLAG} ${DEMVIEWER_LIBS}

riskestimator: $(RISKESTIMATOR_OBJS)
	$(CXX) $^ -o $@ ${OPTFLAG} ${RISKESTIMATOR_LIBS}

# dependencies, blanket
$(OBJS):$(HDRS)


## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(TARGETS)

## end
