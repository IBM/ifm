## Makefile

.PHONY: clean realclean

### these flags can be modified
#OPTFLAG =-O6 -DNDEBUG
#OPTFLAG = -O6
#OPTFLAG =-O2 -funswitch-loops -fpredictive-commoning -fgcse-after-reload -ftree-vectorize
OPTFLAG = -g -O3 t  
OPTFLAG = -g -O3 -qstrict -qsmp=omp 
OPTFLAG = -g -O3 -fopenmp 
#OPTFLAG = -g -qsmp=omp -O3  
#OPTFLAG = -O3 -fno-inline 

### options to change the behavior of the code
#OPTIONS = -DLIM_METHOD1 -DLIM_METHOD2 -DPRE_FILTER -DENABLE_MPI  -DRANGE -D
OPTIONS = -DLIM_METHOD1 -DLIM_METHOD2 -DPRE_FILTER -DPARA
#OPTIONS = -DLIM_METHOD1 -DLIM_METHOD2 -DPRE_FILTER -DENABLE_MPI -DDEBUG 
#OPTIONS = -DLIM_METHOD1 -DLIM_METHOD2 -DPRE_FILTER 
#OPTIONS = -DLIM_METHOD1 -DLIM_METHOD2 -DROUTING_ONLY

### if everything goes fine, should not touch anything below

ifeq ("$(shell uname -m)","i686")
  BITS = -m32
else
  BITS = -m64
endif

CXX           = g++
HDRS          = $(wildcard *.h Core/*.h FileIO/*.h)
OBJS          = $(patsubst %.C,%.o, $(wildcard *.C Core/*.C FileIO/*.C))
NETCDF_LIBS   = -lnetcdf_c++ -lnetcdf

## give the user a chance to override variables
ifneq (,$(wildcard BuildRules/$(USER)@$(HOSTNAME).mk))
  include BuildRules/$(USER)@$(HOSTNAME).mk
endif

CXXFLAGS      = -I. -ICore -IFileIO ${OPTFLAG} ${OPTIONS} ${BITS} ${USER_CXXFLAGS}
LINKFLAG      = -lm -lz ${NETCDF_LIBS} ${USER_LDFLAGS}

## change this line since we have different top level drivers
TARGET = ifm

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $^ -o $@ ${OPTFLAG} ${LINKFLAG}

# dependencies, blanket
$(OBJS):$(HDRS)


## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(TARGET)

## end
