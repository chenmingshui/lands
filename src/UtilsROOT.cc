#include "UtilsROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
void FillTree(TString sfile, vector<int> array){
	TFile fTrees("tree_"+sfile+".root", "RECREATE");
	Double_t brT;
	TTree *tree = new TTree("T","T"); 
	tree->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array.size(); i++){
		brT= array.at(i);
		tree->Fill();
	}
	fTrees.Write();
	fTrees.Close();
	//	if(tree) delete tree;
	//	cout<<"delete me here 6"<<endl;
}
void FillTree(TString sfile, vector<double> array){
	TFile fTrees("tree_"+sfile+".root", "RECREATE");
	Double_t brT;
	TTree *tree = new TTree("T","T"); 
	tree->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array.size(); i++){
		brT= array.at(i);
		tree->Fill();
	}
	fTrees.Write();
	fTrees.Close();
	//	if(tree) delete tree;
	//	cout<<"delete me here 6"<<endl;
}
void FillTree(TString sfile, double * array, int nsize){
	TFile fTrees("tree_"+sfile+".root", "RECREATE");
	Double_t brT;
	TTree *tree = new TTree("T","T"); 
	tree->Branch("brT", &brT, "brT/D");
	for(int i=0; i<nsize; i++){
		brT= array[i];
		tree->Fill();
	}
	fTrees.Write();
	fTrees.Close();
	//	if(tree) delete tree;
	//	cout<<"delete me here 6"<<endl;
}
void FillTree(TString sfile, int* array, int nsize){
	TFile fTrees("tree_"+sfile+".root", "RECREATE");
	Double_t brT;
	TTree *tree = new TTree("T","T"); 
	tree->Branch("brT", &brT, "brT/D");
	for(int i=0; i<nsize; i++){
		brT= array[i];
		tree->Fill();
	}
	fTrees.Write();
	fTrees.Close();
	//	if(tree) delete tree;
	//	cout<<"delete me here 6"<<endl;
}
void FillSigmas(TString sfile, double m2s, double m1s, double mean, double p1s, double p2s, double asimov){
	TFile f("sigmas_"+sfile+".root", "RECREATE");
	TH1F h("h","h", 6, 0.5, 6.5);
	h.SetBinContent(1, m2s);		
	h.SetBinContent(2, m1s);		
	h.SetBinContent(3, mean);		
	h.SetBinContent(4, p1s);		
	h.SetBinContent(5, p2s);		
	h.SetBinContent(6, asimov);		
	h.Write();
	f.Write();
	f.Close();
}

void StringStrip( std::string & str ) {
	//-------------------------------------------------------------------------------
	// Strip spaces at the beginning and the end from a string.
	size_t sPos = 0;
	size_t ePos = str.length();
	while ( str[sPos] == ' ' ) { ++sPos; }
	while ( str[ePos] == ' ' ) { --ePos; }
	str = str.substr( sPos, ePos - sPos );
}

void StringSplit( std::vector < std::string > & splitValues, 
		const std::string & str,
		const std::string & delim ) {
	//-------------------------------------------------------------------------------
	// Split a string by a delimiter and return it's vector of strings.
	std::string str2 = str;

	size_t pos = 0;

	while (( pos = str2.find_first_of( delim )) != std::string::npos ) {
		std::string s = str2.substr(0, pos);
		StringStrip( s );
		if(s!="")splitValues.push_back( s );
		str2 = str2.substr( pos + delim.length());
	}

	StringStrip( str2 );
	if(str2=="")return;
	splitValues.push_back( str2 );
}

TString ReadFile(const char*fileName){
	bool debug = false;

	// open the config file and go through it
	std::ifstream ifile(fileName);

	if(ifile.fail()){
		TString error("File ");
		error+=fileName;
		error+=" could not be opened.";
		cout<<"fReadFile: "<< error.Data()<<endl;
		exit(0);
		//return -1;
	}

	TString ifileContent("");
	ifileContent.ReadFile(ifile);
	ifile.close();

	// Tokenise the file using the "\n" char and parse it line by line to strip
	// the comments.
	TString ifileContentStripped("");

	TObjArray* lines_array = ifileContent.Tokenize("\n");
	TIterator* lineIt=lines_array->MakeIterator();

	bool in_comment=false;
	TString line;
	TObject* line_o;

	while((line_o=(*lineIt)())){ // Start iteration on lines array
		line = (static_cast<TObjString*>(line_o))->GetString();

		// Are we in a multiline comment?
		if (in_comment)
			if (line.EndsWith("*/")){
				in_comment=false;
				//if (fVerbose) Info("fReadFile","Out of multiline comment ...");
				if (debug) cout<<"fReadFile: Out of multiline comment ..."<<endl;

				continue;
			}

		// Was line a single line comment?

		if ((line.BeginsWith("/*") && line.EndsWith("*/")) ||
				line.BeginsWith("//")){
			//if (fVerbose) Info("fReadFile","In single line comment ...");
			if (debug) cout<<"fReadFile: In single line comment ..."<<endl;
			continue;
		}

		// Did a multiline comment just begin?
		if (line.BeginsWith("/*")){
			in_comment=true;
			//if (fVerbose) Info("fReadFile","In multiline comment ...");
			if (debug) printf("fReadFile: In multiline comment ...\n");
			continue;
		}

		ifileContentStripped+=line+"\n";
	}

	delete lines_array;
	delete lineIt;
	in_comment=false;
	return ifileContentStripped;
}
bool ConfigureModel(CountingModel *cms, TString ifileContentStripped){
	bool debug = false;
	// Now proceed with the parsing of the stripped file

	//lines_array = ifileContentStripped.Tokenize(";");
	TObjArray *lines_array = ifileContentStripped.Tokenize("\n");
	TIterator* lineIt=lines_array->MakeIterator();
	vector<TString> lines;
	TObject* line_o;

	const int nNeutrals=2;
	TString neutrals[nNeutrals]={"\t"," "};

	while((line_o=(*lineIt)())){

		TString line = (static_cast<TObjString*>(line_o))->GetString();

		// Strip spaces at the beginning and the end of the line
		line.Strip(TString::kBoth,' ');

		// Put the single statement in one single line
		line.ReplaceAll("\n","");

		// Do we have an echo statement? "A la RooFit"
		if (line.BeginsWith("echo")){
			line = line(5,line.Length()-1);
			if (debug)
				std::cout << "Echoing line " << line.Data() << std::endl;
			std::cout << "[ ] echo: "
				<< line.Data() << std::endl;
			continue;
		}

		if (debug) cout<<"fReadFile: Reading --> "<<line.Data()<<" <--"<<endl;

		if (line == ""){
			if (debug) printf("fReadFile: Empty line: skipping ...\n");
			continue;
		}

		lines.push_back(line);
	}

	delete lineIt;
	delete lines_array;

	//   
	//  For Andrey's input format 
	
	// get number of channels
	int nchannel = 0;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("imax")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(debug)cout<<"NChannel = "<<ss[1]<<endl;
			TString s = ss[1];
			nchannel = s.Atoi();
		} 
	}	
	if(nchannel==0) {
		cout<<"number of channels = 0,  please check imax"<<endl;
	       	return false;
	}

	// get numbers of processes in each channel
	int *nprocesses = new int[nchannel];
	for(int i=0; i<nchannel; i++) nprocesses[i]=0;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("bin")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			for(int i=1; i<ss.size(); i++){
				int bin = (TString(ss[i])).Atoi();
				if(bin<(nchannel+1) && bin>=1) nprocesses[bin-1]++;
				else {
					cout<<"there is a bin number = "<<bin<<", out of [1,"<<nchannel<<"]"<<endl;
					return false;
				}
			}
		} 
	}	
	if(debug){
		for(int i = 0; i<nchannel; i++) cout<<"processes in Channel "<<i<<" = "<<nprocesses[i]<<endl;
	}
	for(int i=0; i<nchannel; i++) if(nprocesses[i]==0) {cout<<"channel "<<i<<", number of processes = 0"<<endl; return false;}

	int ntotprocesses = 0;
	for(int i = 0; i<nchannel; i++) ntotprocesses+=nprocesses[i];

	int *binnumber = new int[ntotprocesses];
	int *subprocess=new int[ntotprocesses];
	int index =0 ;
	for(int c=0; c<nchannel; c++){
		for(int p=0; p<nprocesses[c]; p++) 
		{
			binnumber[index] = (c+1);
			subprocess[index]=p;
			index++;
		}
	}
	if(debug){
		cout<<"bin ";
	       	for(int i=0; i<ntotprocesses; i++) cout<<binnumber[i]<<" ";
		cout<<endl;
	}


	// get expected event rate
	//
	double *eventrate =new double[ntotprocesses];
	for(int i=0; i<ntotprocesses; i++) eventrate[i]=0;
	bool hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("rate")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size()-1<ntotprocesses) {cout<<"number of eventrate < "<<ntotprocesses<<endl; return false;} 
			for(int i=1; i<(ntotprocesses+1); i++){
				float ev = (TString(ss[i])).Atof();
				eventrate[i-1] = ev;
			}
			hasFilled =  true;
		}
	}
	if(hasFilled==false) {cout<<"need a line starting with \"rate\" "<<endl; return false;}
	if(debug){
		cout<<"event rate: ";
		for(int i=0; i<ntotprocesses; i++) cout<<eventrate[i]<<" ";
		cout<<endl;
	}

	// get observed dataset
	//
	double *observeddata=new double[nchannel];
	for(int i=0; i<nchannel; i++) observeddata[i]=0;
	hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("Observation")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size()-1<nchannel) {cout<<"channels of observed data < "<<nchannel<<endl; return false;} 
			for(int i=1; i<(nchannel+1); i++){
				double ev = (TString(ss[i])).Atof();
				observeddata[i-1] = ev;
			}
			hasFilled =  true;
		}
	}
	if(hasFilled==false) {cout<<"need a line starting with \"Observation\" "<<endl; return false;}
	if(debug){
		cout<<"observed data: ";
		for(int i=0; i<nchannel; i++) cout<<observeddata[i]<<" ";
		cout<<endl;
	}

	//  Fill the model with the yields
	index=0;
	for(int c=0; c<nchannel; c++){
		vector<double> sigbkgs; sigbkgs.clear();
		for(int p=0; p<nprocesses[c]; p++) {
			sigbkgs.push_back(eventrate[index]);
			index++;
		}
		cms->AddChannel( sigbkgs );
		cms->AddObservedData(c, observeddata[c]);
	}	

	//cms->Print();

	if(debug) cout<<"filled yields"<<endl;

	//  get number of independant systematics sources
	int nsyssources = 0;
	hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("kmax")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			nsyssources = (TString(ss[1])).Atoi();
			hasFilled = true;
		}
	}
	if (hasFilled == false) cout<<"need kmax to set number of independant systematic sources "<<endl;
	if (nsyssources==0) cout<<"no systematics input"<<endl;

	if(debug) cout<<"number of independant systematics sources = "<<nsyssources<<endl;
	// fill the model with systemics
	for(int s=0; s<nsyssources; s++){
		vector<string>	ss; 
		ss.clear();
		StringSplit(ss, lines[lines.size()-1-s].Data(), " ");
		int indexcorrelation = (TString(ss[0])).Atoi();
		int pdf = 0; 
		if(ss[1]=="lnN") pdf=1;
		else if(ss[1]=="trG") pdf=2;
		else pdf =  (TString(ss[1])).Atoi();

		for(int p=0; p<ntotprocesses; p++){
			double err = (TString(ss[p+2])).Atof()-1.0;
			if(err == 0.) continue;
			//cout<<"delete me: c="<<binnumber[p]-1<<" s="<< subprocess[p]<<endl;
			cms->AddUncertainty(binnumber[p]-1, subprocess[p], err, pdf, indexcorrelation );
		}
	}

	if(debug) cout<<"filled systematics"<<endl;

	return true;
}

bool ConfigureModel(CountingModel *cms, const char* fileName){
	TString s = ReadFile(fileName);
	return ConfigureModel(cms, s);
}

