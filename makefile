LIB  := lands.so
EXEC :=	test/CLs.exe \
	test/Bayesian.exe \
	test/lands.exe

SOURCES := $(wildcard src/*.cc)
OBJECTS := $(patsubst src/%.cc, bin/%.o, $(SOURCES))

ifdef ROOFITSYS
    RF_CFLAGS      := -I${ROOFITSYS}/include
    RF_LINKERFLAGS := -L${ROOFITSYS}/lib
endif

CXX      := g++
CXXFLAGS := -O2 -g -fPIC -Iinclude $(shell root-config --cflags) ${RF_CFLAGS}

LD       := g++
LDFLAGS  := $(shell root-config --libs --ldflags) -lMathMore -lMinuit -lRooFit -lRooFitCore -lFoam ${RF_LINKERFLAGS}

ifeq ($(shell root-config --platform),macosx)
	MACOSXFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -Xlinker -bind_at_load -flat_namespace
endif

# Enable CPU profiling of... requires google performance tools.
# You also need to set CPUPROFILE enviroanment variable.
#PROFILE=1
ifdef PROFILE
  CXXFLAGS += -DLANDS_PROFILE
  LDFLAGS  += -lprofiler
endif


 # stop removing intermediate files
.SECONDARY:

all: $(LIB) $(EXEC)

bin/%.o: src/%.cc
	$(CXX) $(CXXFLAGS)  $< -c -o $@

${LIB}: ${OBJECTS}
	$(LD) $(LDFLAGS) $(MACOSXFLAGS) -shared  ${OBJECTS} -o $@

%.exe: %.cc ${OBJECTS}
	$(LD) $(CXXFLAGS) $(LDFLAGS) -o $@ $< ${OBJECTS}

#for MSSMA
drawMSSMA: test/drawMSSMA.exe

clean:
	@rm -v -f \
	bin/*.o \
	test/*.exe \
	*.so

cleanall: clean
	@rm -v -f \
	test/log* \
	test/*.gif \
	test/*.eps \
	test/*.root \
	*/*~ *~
