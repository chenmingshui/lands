LIB = lands.so
EXEC = test/Significance_dataCard.exe \
	test/MultipleChannels.exe \
	test/CLs.exe \
	test/Bayesian.exe \
	test/ShapeAnalysis.exe \
	test/CLs_dataCard.exe \
	test/Bayesian_dataCard.exe \
	test/Bayesian_2dataCards.exe \
	test/ProfileLikelihoodApproxLimit.exe \
	test/lands.exe

SOURCES = $(wildcard src/*.cc)
COMPONENTS = $(patsubst src%.cc,bin%.o,$(SOURCES))

CC = g++
CFLAGS = -fPIC $(shell root-config --cflags)  -I ./include -I ${ROOTSYS}/include 

LINKER = g++
LINKERFLAGS = $(shell root-config --libs --ldflags) -lMinuit 

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
