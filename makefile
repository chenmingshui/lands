COMPONENTS = bin/BayesianLimitBase.o bin/CRandom.o bin/PdfRandom.o bin/PlotUtilities.o bin/BayesianLimit.o bin/UnbinnedBayesianLimit.o bin/CLsLimit.o bin/Utilities.o bin/UnbinnedCLsLimit.o bin/CountingModel.o bin/UtilsROOT.o bin/BinnedInterface.o  bin/BayesianBase.o bin/LimitBands.o

CFLAGS = -fPIC -I ./include -c -o
CFLAGS1 = -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -c -o

MultipleChannels = test/MultipleChannels.exe
AKinputTemplate = test/AKinputTemplate.exe
CLs = test/CLs.exe
Bayesian = test/Bayesian.exe
ShapeAnalysis = test/ShapeAnalysis.exe

#for MSSMA
drawMSSMA: test/drawMSSMA.exe

AK: test/AKinputTemplate.exe

all: ${CLs} ${MultipleChannels} ${ShapeAnalysis} ${Bayesian} 

test/AKinputTemplate.exe: test/AKinputTemplate.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

test/drawMSSMA.exe: test/drawMSSMA.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

test/ShapeAnalysis.exe: test/ShapeAnalysis.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

test/Bayesian.exe: test/Bayesian.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

test/CLs.exe: test/CLs.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

test/MultipleChannels.exe: test/MultipleChannels.cc ${COMPONENTS}
	g++ -fPIC -I ./include -I /${ROOTSYS}/include `root-config --cflags --libs` -o $@ $< ${COMPONENTS}

bin/BinnedInterface.o: src/BinnedInterface.cc bin/CountingModel.o
	g++ ${CFLAGS1} $@ $<

bin/UtilsROOT.o: src/UtilsROOT.cc bin/CountingModel.o
	g++ ${CFLAGS1} $@ $<

bin/LimitBands.o: src/LimitBands.cc bin/CRandom.o bin/Utilities.o bin/PdfRandom.o bin/CountingModel.o bin/CLsLimit.o bin/BayesianBase.o
	g++ ${CFLAGS} $@ $<

bin/PlotUtilities.o: src/PlotUtilities.cc bin/Utilities.o bin/BayesianLimitBase.o
	g++ ${CFLAGS1} $@ $<

bin/UnbinnedCLsLimit.o: src/UnbinnedCLsLimit.cc bin/CLsLimit.o bin/Utilities.o bin/CRandom.o bin/PdfRandom.o
	g++ ${CFLAGS1} $@ $<

bin/BayesianLimit.o: src/BayesianLimit.cc bin/BayesianLimitBase.o bin/CRandom.o bin/Utilities.o bin/CountingModel.o
	g++ ${CFLAGS} $@ $<

bin/UnbinnedBayesianLimit.o: src/UnbinnedBayesianLimit.cc bin/BayesianLimitBase.o bin/Utilities.o
	g++ ${CFLAGS} $@ $<

bin/CLsLimit.o: src/CLsLimit.cc bin/CRandom.o bin/Utilities.o bin/PdfRandom.o bin/CountingModel.o
	g++ ${CFLAGS} $@ $<

bin/BayesianBase.o: src/BayesianBase.cc bin/Utilities.o bin/CountingModel.o
	g++ ${CFLAGS} $@ $<

bin/BayesianLimitBase.o: src/BayesianLimitBase.cc bin/Utilities.o bin/CountingModel.o
	g++ ${CFLAGS} $@ $<

bin/CountingModel.o: src/CountingModel.cc bin/CRandom.o bin/Utilities.o bin/PdfRandom.o
	g++ ${CFLAGS} $@ $<
	
bin/PdfRandom.o: src/PdfRandom.cc bin/CRandom.o
	g++ ${CFLAGS} $@ $<
	
bin/CRandom.o: src/CRandom.cc
	g++ ${CFLAGS} $@ $<

bin/Utilities.o: src/Utilities.cc
	g++ ${CFLAGS} $@ $<

clean:
	rm bin/*.o  -rf
	rm test/*.exe -rf

cleanall:
	rm bin/*.o -rf
	rm test/*.exe -rf
	rm test/log* -rf
	rm test/*.gif -rf
	rm test/*.eps -rf 
	rm test/*.root -rf
	rm test/*~ -rf 
	rm test/plots -rf
	rm test/roots -rf
	rm test/pilot -rf
	rm */*~ *~ -rf
	
