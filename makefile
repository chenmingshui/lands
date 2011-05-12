LIB = lands.so
EXEC = 	test/CLs.exe \
	test/Bayesian.exe \
	test/lands.exe

SOURCES = $(wildcard src/*.cc)
COMPONENTS = $(patsubst src%.cc,bin%.o,$(SOURCES))

ifdef ROOFITSYS
    RF_CFLAGS = -I ${ROOFITSYS}/include
    RF_LINKERFLAGS = -L ${ROOFITSYS}/lib
endif



CC = g++
CFLAGS = -fPIC $(shell root-config --cflags)  -I ./include -I ${ROOTSYS}/include ${RF_CFLAGS}

LINKER = g++
LINKERFLAGS = $(shell root-config --libs --ldflags) -lMinuit -lRooFit -lRooFitCore -lFoam ${RF_LINKERFLAGS}


ifeq ($(shell root-config --platform),macosx)
	MACOSXFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace
endif

 # stop removing intermediate files
.SECONDARY:

bin/%.o: src/%.cc
	$(CC) $(CFLAGS)  $< -c -o $@
%.so: ${COMPONENTS}
	$(LINKER) $(LINKERFLAGS) $(MACOSXFLAGS) -shared  $(COMPONENTS) -o $@

all: $(LIB) $(EXEC)

%.exe: %.cc $(COMPONENTS)
	$(LINKER) $(CFLAGS) $(LINKERFLAGS) -o $@ $< ${COMPONENTS}

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
