#include "UtilsROOT.h"
#include <memory>
#include <math.h>
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TDirectory.h"
#include "TKey.h"
#include "RooAbsData.h"
map<TString, RooWorkspace*> MAP_RWSname_Pointer;
void SaveResults(TString sfile, double mH, double limit, double limitErr, double significance, double pvalue, double rm2s, double rm1s, double rmedian, double rmean, double rp1s, double rp2s){
	TFile fTrees(sfile+".root", "RECREATE");
	TTree *tree = new TTree("T","T"); 
	tree->Branch("mH", &mH, "mH/D");
	tree->Branch("limit", &limit, "limit/D");
	tree->Branch("limitErr", &limitErr, "limitErr/D");
	tree->Branch("significance", &significance, "significance/D");
	tree->Branch("pvalue", &pvalue, "pvalue/D");
	tree->Branch("rm2s", &rm2s, "rm2s/D");
	tree->Branch("rm1s", &rm1s, "rm1s/D");
	tree->Branch("rmedian", &rmedian, "rmedian/D");
	tree->Branch("rmean", &rmean, "rmean/D");
	tree->Branch("rp1s", &rp1s, "rp1s/D");
	tree->Branch("rp2s", &rp2s, "rp2s/D");
	tree->Fill();

	fTrees.Write();
	fTrees.Close();
}
TH1F* GetHisto(string filename, string histoname){

	//cout<<filename<<", "<<histoname<<endl;
	// FIXME need to check if filename is exist, and histoname is exist 
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; exit(0);};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
	if(f==NULL) f=new TFile(filename.c_str());
	TH1F *h = (TH1F*)f->Get(histoname.c_str());
	//if(!h) {cout<<"hist ["<<histoname<<"] in file ["<<filename<<"] couldn't be found"<<endl; exit(0);};
	return h;
}
TTree * LoadTreeBonly(TString filename, TString & treeName){
	TTree * tb;
	// dedicated to locate a tree with name including Q_data value
	if(treeName!="")return (TTree*)GetTObject(filename.Data(), treeName.Data());
	else{
		if( gSystem->AccessPathName(filename)) {cout<<filename<<" couldn't be found"<<endl; exit(1);};
		TFile *f;
		f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
		if(f==NULL) f=new TFile(filename);
		if (!f) {cout<<"Input file: "<<f->GetName()<<" error!  either not exist or empty,  exit"<<endl; exit(1);}
		TDirectory *toyDir = f->GetDirectory("");
		TIter next(toyDir->GetListOfKeys()); TKey *k;
		int ntrees = 0;
		while ((k = (TKey *) next()) != 0) {
			TString name(k->GetName());
			if(name.BeginsWith("T")){
				name.ReplaceAll("T","");
				cout<<"File ["<<filename<<"] contains tree with -2lnQ for data = "<<name<<endl;
				if(name.IsFloat()==false) { cout<< "     *not float number*  "<<endl ; continue;}
				treeName = name;
				tb = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
				ntrees ++ ;
			}
		}

		if(ntrees>1) cout<<"WARNING: this File contains "<<ntrees<<" trees, we take "<<tb->GetName()<<endl;
		if(ntrees==0) {
			cout<<"ERROR: File ["<<filename<<"] contains 0 tree"<<endl;
			exit(1);
		}
	}

	return tb;

}
TObject* GetTObjectA(string filename, string objname){ // return TObject no matter what you find or not 
	TObject *h;
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; return h;};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
	if(f==NULL) f=new TFile(filename.c_str());
	h = (TObject*)f->Get(objname.c_str());
	//if(!h) {cout<<"object ["<<objname<<"] in file ["<<filename<<"] couldn't be found"<<endl; exit(0);};
	return h;
}
TObject* GetTObject(string filename, string objname){

	cout<<" GetTObject "<<filename.c_str()<<" "<<objname.c_str()<<endl;

	// FIXME need to check if filename is exist, and histoname is exist 
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; exit(0);};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
	if(f==NULL) f=new TFile(filename.c_str());
	TObject *h = (TObject*)f->Get(objname.c_str());
	if(!h) {cout<<"object ["<<objname<<"] in file ["<<filename<<"] couldn't be found"<<endl; exit(0);};
	return h;
}
void FillTree(TString sfile, vector<int> array){
	TFile fTrees(sfile+"_tree.root", "RECREATE");
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
void FillTree(TString sfile, double d1, double d2, vector<double> array1,  vector<double> array2, TString d1Name, TString d2Name, TString array1Name, TString array2Name, TString option ){
	TFile fTrees(sfile, option);
	Double_t brT;
	TTree *tree1 = new TTree(array1Name, array1Name); 
	tree1->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array1.size(); i++){
		brT= array1.at(i);
		tree1->Fill();
	}
	TTree *tree2 = new TTree(array2Name, array2Name); 
	tree2->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array2.size(); i++){
		brT= array2.at(i);
		tree2->Fill();
	}
	TH1D *h1 = new TH1D(d1Name, d1Name, 1, 0, 1);
	h1->SetBinContent(1,d1);
	TH1D *h2 = new TH1D(d2Name, d2Name, 1, 0, 1);
	h2->SetBinContent(1,d2);
	fTrees.Write();
	fTrees.Close();
}
void FillTree2(TString sfile, double d1, double d2, vector<double> array1,  vector<double> array2, TString d1Name, TString array1Name, TString array2Name, TString option ){
	TFile fTrees(sfile, option);
	Double_t brT;
	TTree *tree1 = new TTree(array1Name, array1Name); 
	tree1->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array1.size(); i++){
		brT= array1.at(i);
		tree1->Fill();
	}
	TTree *tree2 = new TTree(array2Name, array2Name); 
	tree2->Branch("brT", &brT, "brT/D");
	for(int i=0; i<array2.size(); i++){
		brT= array2.at(i);
		tree2->Fill();
	}

	TTree *tree3 = new TTree(d1Name, d1Name); 
	tree3->Branch("brT", &brT, "brT/D");
	brT= d2;
	tree3->Fill();
	fTrees.Write();
	fTrees.Close();
}
void FillTree(TString sfile, vector<double> array, TString treeName){
	TFile fTrees(sfile+"_tree.root", "RECREATE");
	Double_t brT;
	TString stmp = "T"; stmp+=treeName;
	TTree *tree = new TTree(stmp, stmp); 
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
	TFile fTrees(sfile+"_tree.root", "RECREATE");
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
	TFile fTrees(sfile+"_tree.root", "RECREATE");
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

bool isWordInMap(TString s, std::map<TString, vector<TString> >tMap){
	std::map<TString, vector<TString> >::iterator p;
	bool bFound=false;
	// Show key
	for(p = tMap.begin(); p!=tMap.end(); ++p)
	{
		TString strKey = p->first;
		if( s == strKey)
		{
			bFound = true;
		}
	}
	return bFound;
}

void StringStrip( std::string & str ) {
	//-------------------------------------------------------------------------------
	// Strip spaces and tab at the beginning and the end from a string.
	size_t sPos = 0;
	size_t ePos = str.length();
	while ( str[sPos] == ' ' || str[sPos]=='\t' ) { ++sPos; }
	while ( str[ePos] == ' ' || str[ePos]=='\t' ) { --ePos; }
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

		string stmp = line.Data();
		StringStrip(stmp);

		line = stmp;

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
				line.BeginsWith("//") || line.BeginsWith("#")){
			//if (fVerbose) Info("fReadFile","In single line comment ...");
			if (debug) cout<<"fReadFile: In single line comment ..."<<endl;
			continue;
		}
		bool bgoodline = false;
		if(line.BeginsWith("-") 
				or line.BeginsWith("0")
				or line.BeginsWith("1")
				or line.BeginsWith("2")
				or line.BeginsWith("3")
				or line.BeginsWith("4")
				or line.BeginsWith("5")
				or line.BeginsWith("6")
				or line.BeginsWith("7")
				or line.BeginsWith("8")
				or line.BeginsWith("9")
				or line.BeginsWith("a")
				or line.BeginsWith("b")
				or line.BeginsWith("c")
				or line.BeginsWith("d")
				or line.BeginsWith("e")
				or line.BeginsWith("f")
				or line.BeginsWith("g")
				or line.BeginsWith("h")
				or line.BeginsWith("i")
				or line.BeginsWith("j")
				or line.BeginsWith("k")
				or line.BeginsWith("l")
				or line.BeginsWith("m")
				or line.BeginsWith("n")
				or line.BeginsWith("o")
				or line.BeginsWith("p")
				or line.BeginsWith("q")
				or line.BeginsWith("r")
				or line.BeginsWith("s")
				or line.BeginsWith("t")
				or line.BeginsWith("u")
				or line.BeginsWith("v")
				or line.BeginsWith("w")
				or line.BeginsWith("x")
				or line.BeginsWith("y")
				or line.BeginsWith("z")
				or line.BeginsWith("A")
				or line.BeginsWith("B")
				or line.BeginsWith("C")
				or line.BeginsWith("D")
				or line.BeginsWith("E")
				or line.BeginsWith("F")
				or line.BeginsWith("G")
				or line.BeginsWith("H")
				or line.BeginsWith("I")
				or line.BeginsWith("J")
				or line.BeginsWith("K")
				or line.BeginsWith("L")
				or line.BeginsWith("M")
				or line.BeginsWith("N")
				or line.BeginsWith("O")
				or line.BeginsWith("P")
				or line.BeginsWith("Q")
				or line.BeginsWith("R")
				or line.BeginsWith("S")
				or line.BeginsWith("T")
				or line.BeginsWith("U")
				or line.BeginsWith("V")
				or line.BeginsWith("W")
				or line.BeginsWith("X")
				or line.BeginsWith("Y")
				or line.BeginsWith("Z")
				)
			       	bgoodline =  true;
		if(bgoodline==false)continue;

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
int GetTotProc(int c, int p, const vector< vector<string> >&vv_procnames){
	int stop = 0;
	for(int i=0; i<c; i++){
		stop += vv_procnames[i].size();
	}
	return stop+p;
}
double GetRate(string c, const vector<string>& channelnames, string p, const vector< vector<string> >& vv_procnames, const vector<string>&ss1){
	int stop = 0;
	for(int i=0; i<channelnames.size(); i++){
		for(int j=0; j< vv_procnames[i].size(); j++){
			if(channelnames[i]==c && vv_procnames[i][j]==p) break;
			stop ++; 
		}
		if(channelnames[i]==c) break;
	}
	TString s = ss1[stop+1];
	if(s == "-") return -9999;
	else return s.Atof();
}
string GetUncertainy(string c, const vector<string> & channelnames, string p, const vector< vector<string> >&vv_procnames,const vector<string>&ss1){
	int stop = 0;
	for(int i=0; i<channelnames.size(); i++){
		for(int j=0; j< vv_procnames[i].size(); j++){
			if(channelnames[i]==c && vv_procnames[i][j]==p) break;
			stop ++; 
		}
		if(channelnames[i]==c) break;
	}
	//cout<<stop+p+3<<endl;
	if(ss1[1]=="gmN" or ss1[1]=="gmA") return ss1[stop+3];
	else return ss1[stop+2];
}
string GetUncertainy(int c, int p, const vector< vector<string> >& vv_procnames, const vector<string>& ss1){
	int stop = 0;
	for(int i=0; i<c; i++){
		stop += vv_procnames[i].size();
	}
	//cout<<stop+p+3<<endl;
	if(ss1[1]=="gmN" or ss1[1]=="gmA") return ss1[stop+p+3];
	else return ss1[stop+p+2];
}
bool CheckIfDoingShapeAnalysis(CountingModel* cms, double mass, TString ifileContentStripped, bool bUseHist, int debug){
	TString smass;
	if(fmod(mass,1)==0) smass.Form("%.0f",mass);
	else smass.Form("%.1f",mass); 
	//int debug = 0;
	vector<TString> lines;
	lines = SplitIntoLines(ifileContentStripped, debug);

	// check if there is key word "shape", if yes, we need expand the
	// whole file to include all bins of shapes
	bool hasShape = false;
	bool hasParametricShape = false;
	vector<TString> shapeinfo; 
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("shapes ")){
			if(GetWordFromLine(lines[l],3)=="FAKE") continue;
			hasShape = true;
			shapeinfo.push_back(lines[l]);
			if(GetWordFromLine(lines[l], 4).Contains(":")) hasParametricShape =true;
		}
	}
	{
		// get number of channels
		int nchannel = -1;
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("imax ")){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				if(debug)cout<<"NChannel = "<<ss[1]<<endl;
				TString s = ss[1];
				if(s.IsDigit()==false) {
					//cout<<"need a number after imax"<<", currently it's not a number but \'"<<s<<"\'"<<endl;
					//exit(0);
					nchannel = -1;
					cout<<"WARNING: You input of imax (number of channels) is not digit, hence LandS will not check the consistency"<<endl;
				}else {
					nchannel = s.Atoi();
				}
			} 
		}	
		if(nchannel==0) {
			cout<<"number of channels = 0,  please check imax"<<endl;
			return false;
		}

		// get observed dataset and cout how many channels in your model
		//
		vector<double> observeddata;
		bool hasFilled = false;
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("Observation ") or lines[l].BeginsWith("observation ")){
				observeddata.clear();
				if(hasFilled) cout<<"WARNING: You have two lines started with \"observation\", we will use the second line"<<endl;
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				for(int i=1; i<ss.size(); i++){
					if(TString(ss[i]).IsFloat()){
						double ev = (TString(ss[i])).Atof();
						observeddata.push_back(ev);
					}
				}
				hasFilled =  true;
			}
		}
		if(hasFilled==false) {cout<<"ERROR: need a line starting with \"Observation\" "<<endl; exit(0);}
		if(nchannel<0) nchannel = observeddata.size(); // if you don't assign a number to imax
		if(nchannel!= observeddata.size()) {
			cout<<"imax = "<<nchannel<<endl;
			cout<<"observeddata.size = "<<observeddata.size()<<endl;
			cout<<"ERROR: number of channels spcified in \"imax\" line is not consistent with \"observation\" line "<<endl;
			exit(0);
		}
		if(debug){
			cout<<"observed data: ";
			for(int i=0; i<nchannel; i++) cout<<observeddata[i]<<" ";
			cout<<endl;
		}


		// there must be one line of "process", which contains enumeration of processes in each channel
		// you might have another line with "process" which contains process names in each channel
		// since in one of the "process" lines,  for each channel, there must be one and only one process enumerated as "0"
		// we will count how many "0" in that line to determine how many channels actually go into your model
		// double check ....
		//
		// first check how many lines started with "process"
		int nlines_with_process = 0;
		vector< vector<string> > vss_processes; vss_processes.clear();
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("process ")){
				nlines_with_process++;
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				vss_processes.push_back(ss);
			}
		}
		if(nlines_with_process<1 or nlines_with_process>2 ){ cout<<"ERROR: you have "<<nlines_with_process<<" lines started with \"process\".  must be one or two"<<endl; exit(0); }
		// read number of channels from process line
		int tmpn = 0;  int l_proc = 0;
		vector<int> nsigproc; nsigproc.clear(); // record number of signal processes in each channel  (with process number <= 0)
		vector<string> processnames; processnames.clear();
		for(int l=0; l<vss_processes.size(); l++){
			for(int i=1; i<vss_processes[l].size(); i++){ // vss_processes[0] is "process"
				if(!TString(vss_processes[l][i]).IsFloat()) break; // you encounter comments or non-numbers, stop processing 
				if(tmpn==nsigproc.size())nsigproc.push_back(0);
				if(TString(vss_processes[l][i]).Atoi()<=0) nsigproc[nsigproc.size()-1]++;
				if(TString(vss_processes[l][i]).Atoi()==0) tmpn++; // as each channel has "0" standing for one of signal processes
				processnames.push_back(vss_processes[l][i]); // we temporarily assign process name as the "number", if another "process" line found,  we'll make change later

				// check that processes number follow   " -2 -1  0  1  2 " ordering 
				int currproc = TString(vss_processes[l][i]).Atoi();
				if(currproc<=0){
					bool procOrderingIsFine = true;
					if(vss_processes[l].size()<=i+1){procOrderingIsFine = false;} 
					else {
						int nextproc = TString(vss_processes[l][i+1]).Atoi();
						if(nextproc!= currproc+1) procOrderingIsFine = false;
					}
					if(procOrderingIsFine == false) {
						cout<<"ERROR: in \"process\" line, the ordering is wrong,  you have to follow the following ordering in each bin: "<<endl;
						cout<<"       .... -2 -1 0 1 2 ....."<<endl;
						exit(0);
					}
				}
			}
			// debug
			for(int i=0; i<nsigproc.size(); i++){
				if(debug)cout<<" n signal process in channel "<<i<<": "<<nsigproc[i]<<endl;
			}


			if(tmpn>0) {
				l_proc = l; // record which line of "process" lines is for enumeration
				break; // already got the wanted line 
			}
		}
		if(nchannel != tmpn ) { cout<<"ERROR: \"observation\" and \"process\" are not consistent"<<endl; exit(0);}
		// read process names 
		if(vss_processes.size()==2){
			for(int i=1; i<vss_processes[1-l_proc].size(); i++){
				if(i>processnames.size()) break;
				processnames[i-1] = vss_processes[1-l_proc][i];
			}
		}

		vector<string> ss_old_rate; 
		ss_old_rate.clear();
		// get numbers of processes in each channel
		int *nprocesses = new int[nchannel];
		for(int i=0; i<nchannel; i++) nprocesses[i]=0;
		bool hasLineStartWithBin = false;
		vector<string> channelnames; channelnames.clear();
		for(int l=0; l<lines.size(); l++){
			//if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==true){
			//	cout<<"WARNING:  there are two lines beginning with \"bin\" in your card. We will take the first one"<<endl;
			//}
			if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==false ){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				if(ss.size() < 1+2*nchannel) continue; // this "bin" line is too short, will process in below code

				// we now can proceed with names of bins instead of just numbers 1, 2, 3, ... N 
				int bin = 0;
				for(int i=1; i<ss.size(); i++){
					if(ss[i]!=ss[i-1]) { bin++; channelnames.push_back(ss[i]);	}
					if(bin<(nchannel+1) && bin>=1) nprocesses[bin-1]++;
				}
				if(bin != nchannel) {cout<<"ERROR: number of channels got from \"bin\" line is not consistent with \"observation\" line"<<endl; exit(0);};

				hasLineStartWithBin = true;
			} 
		}	
		if(!hasLineStartWithBin) {
			cout<<"Line beginning with \"bin\" is not found in your card"<<endl;
			exit(0);
		}
		if(debug){
			for(int i = 0; i<nchannel; i++) cout<<"processes in Channel "<<i<<" = "<<nprocesses[i]<<endl;
		}
		for(int i=0; i<nchannel; i++) if(nprocesses[i]==0) {cout<<"ERROR: channel "<<i<<", number of processes = 0"<<endl; exit(0);}

		// if you have "bins" or "binname" line,  then replace channel name  with the ones from this line
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("bins ") or lines[l].BeginsWith("binname ") or lines[l].BeginsWith("bin ") ){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				if(ss.size() >= 1+nchannel*2) continue;// this "bin" line is too long, processed in above code 
				for(int i=1; i<ss.size(); i++){
					if(i>nchannel) break;
					channelnames[i-1]=ss[i];
				}
			}
		}

		int ntotprocesses = 0;
		for(int i = 0; i<nchannel; i++) ntotprocesses+=nprocesses[i];

		int *binnumber = new int[ntotprocesses];
		int *subprocess=new int[ntotprocesses];
		int index =0 ;
		vector< vector<string> > vv_procnames; vv_procnames.clear();
		for(int c=0; c<nchannel; c++){
			vector<string> vtmp; vtmp.clear();
			for(int p=0; p<nprocesses[c]; p++) 
			{
				binnumber[index] = (c+1);
				subprocess[index]=p;
				vtmp.push_back(processnames[index]);
				index++;
			}
			vv_procnames.push_back(vtmp);
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
		hasFilled = false;
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("rate ")){
				TString tmps = TString::Format("%10s ","rate");
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				if(ss.size()-1<ntotprocesses) {cout<<"number of eventrate < "<<ntotprocesses<<endl; exit(0);} 
				for(int i=1; i<(ntotprocesses+1); i++){
					float ev = (TString(ss[i])).Atof();
					eventrate[i-1] = ev;
					tmps+= TString::Format("%7.2f ",ev);
				}
				hasFilled =  true;
				StringSplit(ss_old_rate, lines[l].Data(), " ");
			}
		}
		if(hasFilled==false) {cout<<"need a line starting with \"rate\" "<<endl; exit(0);}
		if(debug){
			cout<<"event rate: ";
			for(int i=0; i<ntotprocesses; i++) cout<<eventrate[i]<<" ";
			cout<<endl;
		}


		// FIXME we won't need this "kmax" keyword if we don't want to do a sanity check  ....
		//  get number of independant systematics sources
		int kmax = -1;
		//hasFilled = false;
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("kmax ")){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, lines[l].Data(), " ");
				if(ss.size()>=2 && TString(ss[1]).IsFloat()){
					kmax = (TString(ss[1])).Atoi();
				}else{
					cout<<"WARNING: You input of kmax (number of systematic lines ) is not digit, hence LandS will not check the consistency"<<endl;
					kmax = -1;
				}
				//hasFilled = true;
			}
		}
		//if (hasFilled == false) cout<<"need \"kmax\" to set number of independant systematic sources "<<endl;
		if (kmax==0) cout<<"no systematics input, kmax=0"<<endl;


		// fill the model with systemics
		// sections are separated by "---------"  in the data card
		// last section are all uncertainties,  it's analyzers' responsibility to make it right  
		int nsyssources = 0;
		vector<string> uncernames; uncernames.clear();
		vector<string> uncertypes; uncertypes.clear();
		vector< vector<string> >uncerlines; uncerlines.clear();
		vector< vector<string> >uncerlinesAffectingShapes; uncerlinesAffectingShapes.clear();
		for(int j=lines.size()-1; j>=0; j--){
			if(lines[j].BeginsWith("---") )	break;
			if(lines[j].BeginsWith("rate "))	break;
			if(lines[j].BeginsWith("bin ") )	break;
			if(lines[j].BeginsWith("process "))break;
			if(lines[j].BeginsWith("imax "))break;
			if(lines[j].BeginsWith("jmax "))break;
			if(lines[j].BeginsWith("kmax "))break;
			if(lines[j].BeginsWith("shapes "))continue;
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[j].Data(), " ");
			if( ss.size()<ntotprocesses+2 and (GetWordFromLine(lines[j], 1)!="param" and GetWordFromLine(lines[j],1)!="flatParam") ){
				cout<<"uncertainty configuration is not correct"<<endl; 
				cout<<lines[j]<<endl;
				exit(0);
			}else if(ss[1]=="affects"){
				uncerlinesAffectingShapes.push_back(ss);
			}else{
				nsyssources++;
				uncernames.push_back(ss[0]);
				uncertypes.push_back(ss[1]);
				uncerlines.push_back(ss);
			}
		}
		if(kmax>=0 && kmax!=nsyssources) {cout<<"kmax !=  number of independant uncertainties"<<endl; exit(0);}

		if(debug) cout<<"number of independant systematics sources = "<<nsyssources<<endl;

		if(debug) cout<<"*************shape analysis**************"<<endl;
		TString cardExpanded;
		vector<TString> newlines; // will modify old file line by line

		// must eliminate the "shape" in this new card,
		// otherwise it will enter into unlimited loop

		if(debug) cout<<"shapeinfo lines = "<<shapeinfo.size()<<endl;

		int INDEXofProcess = 1, INDEXofChannel=2;
		vector< vector<string> > shape;
		vector< vector<string> > shapeuncertainties;
		vector< vector<string> > xx_needtointerpret;
		vector< vector<string> > parametricShapeLines; 
		vector< vector<string> > histShapeLines; 
		for(int j=0; j<shapeinfo.size(); j++){
			if(shapeinfo[j].BeginsWith("shapes ")){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, shapeinfo[j].Data(), " ");
				if(ss.size()>=5){
					if(ss[1]=="*" || ss[2]=="*") xx_needtointerpret.push_back(ss);
					else {
						string n4 = ss[4];
						string n5 ; if(ss.size()>5) n5= ss[5];
						ss[4] = TString(ss[4]).ReplaceAll("$MASS", smass);
						ss[4] = TString(ss[4]).ReplaceAll("$CHANNEL", ss[2]);
						ss[4] = TString(ss[4]).ReplaceAll("$PROCESS", ss[1]);
						if(ss.size()>5){
						ss[5] = TString(ss[5]).ReplaceAll("$MASS", smass);
						ss[5] = TString(ss[5]).ReplaceAll("$CHANNEL", ss[2]);
						ss[5] = TString(ss[5]).ReplaceAll("$PROCESS", ss[1]);
						}
						shape.push_back(ss);
						ss[4] = n4; if(ss.size()>5) ss[5]=n5;
						if(ss.size()>5){
							if(ss.size()==6){
								string ss6=ss[5];
								ss.push_back(ss6);

								// for channel c, process t,  looking for uncer which is shaping 
								for(int u=0; u<nsyssources; u++){
									TString type = uncerlines[u][1];
									if(type.BeginsWith("shape")){ 
										ss[4]=uncerlines[u][0];
										string n5 = ss[5];
										string n6 = ss[6];

										TString unc = GetUncertainy(ss[INDEXofChannel], channelnames, ss[INDEXofProcess], vv_procnames, uncerlines[u]);
										if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
											ss[5] = TString(ss[5]).ReplaceAll("$SYSTEMATIC", ss[4]+"Up");
											ss[6] = TString(ss[6]).ReplaceAll("$SYSTEMATIC", ss[4]+"Down");
											ss[5] = TString(ss[5]).ReplaceAll("$MASS", smass);
											ss[6] = TString(ss[6]).ReplaceAll("$MASS", smass);
											ss[5] = TString(ss[5]).ReplaceAll("$CHANNEL", ss[2]);
											ss[5] = TString(ss[5]).ReplaceAll("$PROCESS", ss[1]);
											ss[6] = TString(ss[6]).ReplaceAll("$CHANNEL", ss[2]);
											ss[6] = TString(ss[6]).ReplaceAll("$PROCESS", ss[1]);
											shapeuncertainties.push_back(ss);
											ss[5] = n5; ss[6]=n6;
										}
									}
								}
							}else shapeuncertainties.push_back(ss);
						}
					}
				}
			}
		}
		if(debug)cout<<"xx_needtointerpret.size="<<xx_needtointerpret.size()<<endl;
		if(debug)cout<<"shape.size="<<shape.size()<<endl;
		if(debug)cout<<"shapeuncertainties.size="<<shapeuncertainties.size()<<endl;

		for(int i=0; i<xx_needtointerpret.size(); i++){
			if(xx_needtointerpret[i][INDEXofChannel]!="*" or xx_needtointerpret[i][INDEXofProcess]!="*"){
				for(int c=0; c<nchannel; c++){
					if(xx_needtointerpret[i][INDEXofChannel]=="*" or xx_needtointerpret[i][INDEXofChannel]==channelnames[c]){
						for(int p=0; p<vv_procnames[c].size(); p++){
							if(xx_needtointerpret[i][INDEXofProcess]=="*" or xx_needtointerpret[i][INDEXofProcess]==vv_procnames[c][p]){
								vector<string> newline = xx_needtointerpret[i];
								newline[INDEXofChannel] = channelnames[c]; newline[INDEXofProcess]=vv_procnames[c][p];
								newline[4] = TString(newline[4]).ReplaceAll("$CHANNEL", channelnames[c]);
								newline[4] = TString(newline[4]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
								newline[4] = TString(newline[4]).ReplaceAll("$MASS", smass);
								shape.push_back(newline);
								if(newline.size()>5){
									if(debug) cout<<"debug 0"<<endl;
									if(newline.size()==6){
										string ss6=newline[5];
										newline.push_back(ss6);
										// for channel c, process t,  looking for uncer which is shaping 
										for(int u=0; u<nsyssources; u++){
											TString type = uncerlines[u][1];
											if(type.BeginsWith("shape")){ 
												if(debug) cout<<"debug 1"<<endl;
												newline[4]=uncerlines[u][0];
												string n5 = newline[5];
												string n6 = newline[6];

												TString unc = GetUncertainy(c, p, vv_procnames, uncerlines[u]);
												if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
													newline[5] = TString(newline[5]).ReplaceAll("$SYSTEMATIC", newline[4]+"Up");
													newline[6] = TString(newline[6]).ReplaceAll("$SYSTEMATIC", newline[4]+"Down");
													newline[5] = TString(newline[5]).ReplaceAll("$CHANNEL", channelnames[c]);
													newline[5] = TString(newline[5]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
													newline[6] = TString(newline[6]).ReplaceAll("$CHANNEL", channelnames[c]);
													newline[6] = TString(newline[6]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
													newline[5] = TString(newline[5]).ReplaceAll("$MASS", smass);
													newline[6] = TString(newline[6]).ReplaceAll("$MASS", smass);
													shapeuncertainties.push_back(newline);
													newline[5] = n5; newline[6]=n6;
												}
											}
										}
									}
									else shapeuncertainties.push_back(newline);
								}
							}
						}
						if(xx_needtointerpret[i][INDEXofProcess]=="*" or xx_needtointerpret[i][INDEXofProcess]=="data_obs"){
							vector<string> newline = xx_needtointerpret[i];
							newline[INDEXofChannel] = channelnames[c]; newline[INDEXofProcess]="data_obs";
							newline[4] = TString(newline[4]).ReplaceAll("$CHANNEL", channelnames[c]);
							newline[4] = TString(newline[4]).ReplaceAll("$PROCESS", "data_obs");
							newline[4] = TString(newline[4]).ReplaceAll("$MASS", smass);
							shape.push_back(newline);
						}
					}
				}
			}
		}
		for(int i=0; i<xx_needtointerpret.size(); i++){
			if(xx_needtointerpret[i][INDEXofChannel]=="*" and xx_needtointerpret[i][INDEXofProcess]=="*"){
			//if(xx_needtointerpret[i][INDEXofChannel]!="*" or xx_needtointerpret[i][INDEXofProcess]!="*"){
				for(int c=0; c<nchannel; c++){
					if(xx_needtointerpret[i][INDEXofChannel]=="*" or xx_needtointerpret[i][INDEXofChannel]==channelnames[c]){
						for(int p=0; p<vv_procnames[c].size(); p++){
							if(xx_needtointerpret[i][INDEXofProcess]=="*" or xx_needtointerpret[i][INDEXofProcess]==vv_procnames[c][p]){
								vector<string> newline = xx_needtointerpret[i];
								newline[INDEXofChannel] = channelnames[c]; newline[INDEXofProcess]=vv_procnames[c][p];
								newline[4] = TString(newline[4]).ReplaceAll("$CHANNEL", channelnames[c]);
								newline[4] = TString(newline[4]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
								newline[4] = TString(newline[4]).ReplaceAll("$MASS", smass);
								shape.push_back(newline);
								if(newline.size()>5){
									if(debug) cout<<"debug 0"<<endl;
									if(newline.size()==6){
										string ss6=newline[5];
										newline.push_back(ss6);
										// for channel c, process t,  looking for uncer which is shaping 
										for(int u=0; u<nsyssources; u++){
											TString type = uncerlines[u][1];
											if(type.BeginsWith("shape")){ 
												if(debug) cout<<"debug 1"<<endl;
												newline[4]=uncerlines[u][0];
												string n5 = newline[5];
												string n6 = newline[6];

												TString unc = GetUncertainy(c, p, vv_procnames, uncerlines[u]);
												if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
													newline[5] = TString(newline[5]).ReplaceAll("$SYSTEMATIC", newline[4]+"Up");
													newline[6] = TString(newline[6]).ReplaceAll("$SYSTEMATIC", newline[4]+"Down");
													newline[5] = TString(newline[5]).ReplaceAll("$CHANNEL", channelnames[c]);
													newline[5] = TString(newline[5]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
													newline[6] = TString(newline[6]).ReplaceAll("$CHANNEL", channelnames[c]);
													newline[6] = TString(newline[6]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
													newline[5] = TString(newline[5]).ReplaceAll("$MASS", smass);
													newline[6] = TString(newline[6]).ReplaceAll("$MASS", smass);
													shapeuncertainties.push_back(newline);
													newline[5] = n5; newline[6]=n6;
												}
											}
										}
									}
									else shapeuncertainties.push_back(newline);
								}
							}
						}
						if(xx_needtointerpret[i][INDEXofProcess]=="*" or xx_needtointerpret[i][INDEXofProcess]=="data_obs"){
							vector<string> newline = xx_needtointerpret[i];
							newline[INDEXofChannel] = channelnames[c]; newline[INDEXofProcess]="data_obs";
							newline[4] = TString(newline[4]).ReplaceAll("$CHANNEL", channelnames[c]);
							newline[4] = TString(newline[4]).ReplaceAll("$PROCESS", "data_obs");
							newline[4] = TString(newline[4]).ReplaceAll("$MASS", smass);
							shape.push_back(newline);
						}
					}
				}
			}
		}
		// remove duplicate channel/process 
		for(int i=0; i<shape.size(); i++){
			for(int j=i+1; j<shape.size(); j++){
				if(shape[j][1]==shape[i][1] and shape[j][2]==shape[i][2]) {
					shape[j][0] = "##";
					shape[j][1] = "##";
					shape[j][2] = "##";
				}
			}
		}

		for(int i=0; i<shape.size(); i++){
			if(TString(shape[i][4]).Contains(":") and TString(shape[i][0]) != "##") parametricShapeLines.push_back(shape[i]);
			else if(bUseHist && TString(shape[i][0])!="##") histShapeLines.push_back(shape[i]);
		}
		for(int i=0; i<shapeuncertainties.size(); i++){// for histogram morphing ....    keyword:  shapeN, shapeL, shape, shapeQ ,  alternative two sets of histogram for the variations w.r.t a particular source of systematics
//need to implement shape2 and shapeN2  which do interpolation in log scale   instead of  linear scale 
			for(int j=i+1; j<shapeuncertainties.size(); j++){
				TString stmpi = shapeuncertainties[i][5], stmpj = shapeuncertainties[j][5];
				stmpi.ReplaceAll(smass,"");
				stmpj.ReplaceAll(smass,"");
				if(shapeuncertainties[j][1]==shapeuncertainties[i][1] and shapeuncertainties[j][2]==shapeuncertainties[i][2]
						and (shapeuncertainties[j][5] == shapeuncertainties[i][5] or stmpi==stmpj)
				  ) {
					shapeuncertainties[j][0] = "##";
					shapeuncertainties[j][1] = "##";
					shapeuncertainties[j][2] = "##";
				}
			}
		}

		if(debug){
			cout<<"shape.size="<<shape.size()<<endl;
			for(int i=0; i<shape.size(); i++){
				for(int j=0; j<shape[i].size(); j++){
					cout<<shape[i][j]<<" ";
				}
				cout<<endl;
			}
			if(debug)cout<<"uncertainties.size="<<shapeuncertainties.size()<<endl;
			for(int i=0; i<shapeuncertainties.size(); i++){
				for(int j=0; j<shapeuncertainties[i].size(); j++){
					cout<<shapeuncertainties[i][j]<<" ";
				}
				cout<<endl;
			}
		}

		// need to make a lot of checks:
		// same binning within a channel 
		// normalization of histogram should be equal to expected rate/total observation number
		// normalize shape uncertainties to 1. ,  and also record normalization of  shift_up, shift_down

		if(bUseHist){
			shape.clear();
		}

		int newchannels = nchannel;
		for(int c = 0; c<channelnames.size(); c++){
			for(int j=0; j<shape.size(); j++){
				TString s = shape[j][INDEXofChannel];
				if(s!=channelnames[c])continue;
				if(TString(shape[j][4]).Contains(":")) continue; // parametric channel, will process later ////  vhbb has put histogram into RooDataHist in RooWorkspace, so game changes here 
// to much work to implement for them, so workaround is to convert RooDataHist into TH1F and put them into same TFiles ...
				TH1F *h ;
				/*if(TString(shape[j][4]).Contains(":")) {
					RooWorkspace * w = (RooWorkspace*) GetTObject(shape[j][3], GetWordFromLine(lines[i][4], 0, ":").Data());
					RooDataHist * rdh = (RooDataHist*) 
				} else {
					h= (TH1F*)GetHisto(shape[j][3],shape[j][4]);
				}
				*/
				h= (TH1F*)GetHisto(shape[j][3],shape[j][4]);
				if(h==NULL) continue;
				newchannels+=h->GetNbinsX();
				newchannels-=1;
				delete h;
				break;
			}
		}

		if(debug) cout<<"new number of channels: "<<newchannels<<endl;

		TString s = "imax "; s+=newchannels;
		newlines.push_back(s);
		s = "jmax *"; 
		newlines.push_back(s);
		s = "kmax "; s+=nsyssources;
		newlines.push_back(s);

		for(int k=0; k<lines.size(); k++){
			if(lines[k].BeginsWith("Observation ") or lines[k].BeginsWith("observation ")){
				vector<string> ss_old; ss_old.clear();
				StringSplit(ss_old, lines[k].Data(), " ");
				TString s = "Observation ";
				int n = 0; // for index of new whole channels
				for(int c=0; c<nchannel; c++){//old channel number will be replaced 
					bool isShapeChannel = false;
					for(int p=0; p<shape.size(); p++){
						if(shape[p][INDEXofChannel] !=channelnames[c] || shape[p][INDEXofProcess]!="data_obs") continue; // find the signal one
						if(TString(shape[p][4]).Contains(":")) continue; // parametric channel, will process later
						isShapeChannel=true;
						TH1F* h = (TH1F*)GetHisto(shape[p][3], shape[p][4]);
						if(h==NULL) continue;
						for(int r=1; r<=h->GetNbinsX(); r++){
							n++; s+=h->GetBinContent(r); s+=" ";
						}
						delete h;
					}
					if(isShapeChannel==false) {n++; s+=ss_old[c+1]; s+=" ";}
				}	
				newlines.push_back(s);
			}
		}

		if(debug) cout<<"refill shape with observation"<<endl;

		TString s_bin = "bin ";
		TString s_rate = "rate ";
		TString s_process = "process ";
		TString s_process_name = "process ";
		vector<TString> vs_unc; vs_unc.clear();
		for(int u=0; u<nsyssources; u++) {
			vector<string>	ss1 = uncerlines[u]; 
			if(ss1[1]=="shapeN" or ss1[1]=="shapeN2") ss1[1]="lnN"; 
			TString s = ""; s+=ss1[0].c_str(); s+=" "; s+=ss1[1].c_str(); s+=" "; 
			if(ss1[1]=="gmA" or ss1[1]=="gmN"){
				s+=ss1[2]; s+=" ";
				//cout<<s<<endl;
			};
			//if(uncerlines[u][1]=="gmN") cout<<s<<endl;
			if(uncertypes[u]=="param" || uncertypes[u]=="flatParam") {
				for(int ii = 2; ii<ss1.size(); ii++){
					s+=ss1[ii]; s+=" ";
				}
			}
			vs_unc.push_back(s);
		}

		vector< vector<string> > shape_ch, uncer_ch;	
		int n = 0; // for index of new whole channels
		for(int c=0; c<nchannel; c++){//old channel number will be replaced 
			bool isShapeChannel = false;
			for(int p=0; p<shape.size(); p++){
				if(shape[p][INDEXofChannel]!=channelnames[c] || shape[p][INDEXofProcess]!="data_obs") continue; // find the signal one
				if(TString(shape[p][4]).Contains(":")) continue; // parametric channel, will process later

				shape_ch.clear(); uncer_ch.clear();
				for(int q=0; q<shape.size(); q++){
					if(shape[q][INDEXofChannel]==channelnames[c]) shape_ch.push_back(shape[q]);
				}
				for(int q=0; q<shapeuncertainties.size(); q++){
					if(shapeuncertainties[q][INDEXofChannel]==channelnames[c]) uncer_ch.push_back(shapeuncertainties[q]);
				}

				isShapeChannel=true;

				if(debug)cout<<"channel="<<shape[p][INDEXofChannel]<<" data_obs file="<<shape[p][3]<<" hist="<<shape[p][4]<<endl;

				int n_proc = nprocesses[c];
				TH1F* h = (TH1F*)GetHisto(shape[p][3], shape[p][4]);
				if(h==NULL) continue;
				TH1F* hn[n_proc];
				TH1F *hnorm[n_proc]; // normalized to 1 
				TH1F *hunc_up[n_proc][nsyssources];
				TH1F *hunc_dn[n_proc][nsyssources];
				TH1F *hunc_up_norm[n_proc][nsyssources];
				TH1F *hunc_dn_norm[n_proc][nsyssources];
				for(int x=0; x<n_proc; x++) { 
					hn[x]=NULL; hnorm[x]=NULL; 
					for(int y=0;y<nsyssources;y++){
						hunc_up[x][y]=NULL;
						hunc_dn[x][y]=NULL;
						hunc_up_norm[x][y]=NULL;
						hunc_dn_norm[x][y]=NULL;
					}
				}
				for(int t=0; t<n_proc; t++){
					for(int q=0; q<shape_ch.size(); q++){
						if(shape_ch[q][INDEXofProcess]!=vv_procnames[c][t]) continue; 
						TH1F *htmp =((TH1F*)GetHisto(shape_ch[q][3], shape_ch[q][4]));
						if(htmp==NULL) continue;
						hn[t] = htmp;
		
						if(hn[t]->Integral() == 0) { 
							cout<<"WARNING: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
								<<"] is shape, but the norminal histogram->Integral = 0"<<endl;
							//exit(0);
						}

						// check if normalization consistent with the rate declared in the data card
						double tmp_rate = GetRate(channelnames[c], channelnames, vv_procnames[c][t], vv_procnames, ss_old_rate);
						if(tmp_rate<0) {} //don't check normalization of histogram,  users' responsibility 
						else{ 
							if(fabs(tmp_rate-hn[t]->Integral())/hn[t]->Integral() > 1){ // if difference > 100%, then put Error and exit
								cout<<"ERROR: in channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]<<"], ";
								cout<<" declared rate in data card = "<< tmp_rate<<", not consistent with the integral of input histogram,";
							        cout<<" which is "<<hn[t]->Integral()<<endl;
								cout<<" NOTE: if you don't want to check the consistency, you can put \"-\" in the datacard for the rate. "<<endl;
								cout<<"       In that case, you have the responsibility to make sure the integral (normalization) of histogram makes sense."<<endl;
								exit(0);
							}
						}
							

						TString sname  = "del_clone"; sname+=t;
						hnorm[t]=(TH1F*)hn[t]->Clone(sname);
						if(hn[t]->Integral()!=0)hnorm[t]->Scale(1./hn[t]->Integral());
					}
					if(hnorm[t]==NULL) {
						cout<<" channle ["<<channelnames[c]<<"] proc ["<<vv_procnames[c][t]<<"]"<<" histogram is not found "<<endl;
						exit(1);
					}

					// for channel c, process t,  looking for uncer which is shaping 
					for(int u=0; u<nsyssources; u++){
						TString type = uncerlines[u][1];
						if(type.BeginsWith("shape")){ 
							TString unc = GetUncertainy(c, t, vv_procnames, uncerlines[u]);
							if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
								bool filled = false;
								bool bDEBUGGING = false;
								for(int i=0; i<uncer_ch.size(); i++){ // look for the histograms from the list 
									//ch1_ge3t] process [ttH] sys [CMS_res_j
									if(
											uncer_ch[i][INDEXofProcess]==vv_procnames[c][t]
											and uncer_ch[i][4]==uncerlines[u][0]){
										TH1F * htmp = (TH1F*)GetHisto(uncer_ch[i][3], uncer_ch[i][5]);
										if(htmp==NULL) continue;
										hunc_up[t][u]  = htmp ;
										if(hunc_up[t][u]->Integral()== 0 && hnorm[t]->Integral()!=0) { 
											cout<<" hnorm "<<hnorm[t]->GetName()<<": integral = "<<hnorm[t]->Integral()<<endl;
											cout<<"ERROR: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
												<<"] is shape, but the up_shift histogram->Integral = 0"<<endl;
											exit(0);
										}
										htmp = (TH1F*)GetHisto(uncer_ch[i][3], uncer_ch[i][6]);
										if(htmp==NULL) continue;
										hunc_dn[t][u] =  htmp;
										if(hunc_dn[t][u]->Integral()== 0 && hnorm[t]->Integral()!=0) { 
											cout<<"ERROR: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
												<<"] is shape, but the down_shift histogram->Integral = 0"<<endl;
											exit(0);
										}
										TString sname  = "up_clone"; sname+=t; sname+=u;
										hunc_up_norm[t][u]=(TH1F*)hunc_up[t][u]->Clone(sname);
										if(hunc_up[t][u]->Integral()!=0)hunc_up_norm[t][u]->Scale(1./hunc_up[t][u]->Integral());
										sname  = "dn_clone"; sname+=t; sname+=u;
										hunc_dn_norm[t][u]=(TH1F*)hunc_dn[t][u]->Clone(sname);
										if(hunc_dn[t][u]->Integral()!=0)hunc_dn_norm[t][u]->Scale(1./hunc_dn[t][u]->Integral());
										if(hnorm[t]==NULL) {
											cout<<" channle ["<<channelnames[c]<<"] proc ["<<vv_procnames[c][t]<<"]"<<" histogram is not found "<<endl;
											
											exit(1);
										}
										if(hnorm[t]->Integral()==0){
											cout<<" DEBUGGGGGGG 2"<<endl;
											for(int ii=0; ii<=hunc_up_norm[t][u]->GetNbinsX(); ii++){
											hunc_dn_norm[t][u]->SetBinContent(ii,0);
											hunc_up_norm[t][u]->SetBinContent(ii,0);
											}
										}
										if(0&&debug){
											cout<<"channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]<<": sys."<<uncerlines[u][0]<<endl;
											cout<<hnorm[t]->GetName()<<":"<<hnorm[t]->Integral()<<endl;
											cout<<hunc_dn_norm[t][u]->GetName()<<":"<<hunc_dn_norm[t][u]->Integral()<<endl;
											cout<<hunc_up_norm[t][u]->GetName()<<":"<<hunc_up_norm[t][u]->Integral()<<endl;
										}
										filled = true;
									}
								}	
								if(!filled and !type.Contains("?")) {
									cout<<"ERROR channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]<<"] sys ["<<uncerlines[u][0]
										<<"] is shape uncertainty, we need corresponding histogram "<<endl;
									exit(0);
								}
							}
						}
					}
				}

				if(debug) cout<<" processing "<<h->GetName()<<" "<<channelnames[c]<<"  nbin="<<h->GetNbinsX()<<endl;


				for(int r=1; r<=h->GetNbinsX(); r++){
					n++; 
					int proc = 0;
					for(int t=0; t<n_proc; t++){
						//s_bin+=channelnames[c]; s_bin+="_"; s_bin+=r; s_bin+=" ";
						s_bin+="TH1F_"; s_bin+=channelnames[c];
						s_bin+="_"; s_bin+=h->GetBinLowEdge(r);
						s_bin+="_"; s_bin+=h->GetBinCenter(r); 
						s_bin+="_"; s_bin+=h->GetBinWidth(r); 
						s_bin+=" ";
						s_rate+=((hn[proc]))->GetBinContent(r);s_rate+=" ";
						// fill  s_process  //FIXME
						s_process += (proc-nsigproc[c]+1) ; s_process += " ";
						s_process_name += vv_procnames[c][proc] ; s_process_name += " ";

						if(debug>=10)cout<<"channel "<<c <<", histo bin"<<r<<" oldbin "<<t<<endl;

						for(int u=0; u<nsyssources; u++){
							//for(int k=0; k<uncerlines.size(); k++)
							if(uncertypes[u]=="param" or uncertypes[u]=="flatParam") continue;
							if(1){
								if(debug>=100) cout<<"debug "<<uncertypes[u]<<endl;
								if(uncertypes[u]=="shape" or uncertypes[u]=="shape2" or uncertypes[u]=="shapeL" or uncertypes[u]=="shapeQ" or uncertypes[u]=="shapeN" or uncertypes[u]=="shapeN2" or uncertypes[u]=="shapeStat" or (uncertypes[u]=="shape?" and hunc_dn_norm[t][u]!=NULL)){
									TString unc = GetUncertainy(c, t, vv_procnames, uncerlines[u]);
									if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
										if(hunc_dn_norm[t][u]==NULL){
											cout<<" DEBUGGGGG 1 "<<uncertypes[u]<<"   "<<endl;
											cout<<uncerlines[u][0]<<endl;
											cout<<" hnorm "<<hnorm[t]->GetName()<<endl;
											cout<<" "<<channelnames[c]<<" "<<vv_procnames[c][proc]<<endl;
											cout<<" hunc_dn_norm "<<t<<" "<<u<<" exist "<<endl;
										}
										double down = hunc_dn_norm[t][u]->GetBinContent(r);
										if(debug>10) cout<<"down "<<down<<endl;
										double up = hunc_up_norm[t][u]->GetBinContent(r);
										double norminal = hnorm[t]->GetBinContent(r); 
										if(down>-9999999 && down<9999999 && up>-9999999 && up<99999999 && norminal>-9999999 && norminal<99999999) {
											// ok
										}else{
											cout<<" input not ok:  bin content is crazy : "<<endl;
											cout<<"down = "<<down<<endl;
											cout<<"up   = "<<up<<endl;
											cout<<"norminal = "<<norminal<<endl;
											cout<<hunc_dn_norm[t][u]->GetName()<<":"<<hunc_dn_norm[t][u]->Integral()<<endl;
											cout<<hunc_up_norm[t][u]->GetName()<<":"<<hunc_up_norm[t][u]->Integral()<<endl;
											exit(1);
										}
										if(norminal<0){ 
											// make sure no negativae bin content
											cout<<"WARNING: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]<<" has bin with content < 0 : "<<norminal<<".   we set it to be 0"<<endl; 
											norminal=0;
										}
										double norm_down = hunc_dn[t][u]->Integral();
										double norm_up = hunc_up[t][u]->Integral();
										double norm_norminal = hn[t]->Integral();
										TString stmp;
										if(uncertypes[u]=="shapeN" or uncertypes[u]=="shapeN2" or uncertypes[u]=="shapeStat"){
											vs_unc[u]+=(norminal==0?1:down/norminal); vs_unc[u] += "/";
											stmp.Form("%g",(norminal==0?1:up/norminal));
											//vs_unc[u]+=(norminal==0?1:up/norminal); vs_unc[u] += " ";
											vs_unc[u]+=stmp; vs_unc[u] += " ";
											if(debug>=10)cout<<vs_unc[u]<<endl;
										}else {
											vs_unc[u]+=down; vs_unc[u] += "/";
											stmp.Form("%g",up);
											vs_unc[u]+=stmp; vs_unc[u] += "/";
											stmp.Form("%g",norminal);
											vs_unc[u]+=stmp; vs_unc[u] += "/";
											stmp.Form("%g",norm_norminal);
											vs_unc[u]+=stmp; vs_unc[u] += "/";
											stmp.Form("%g",norm_down/norm_norminal);
											vs_unc[u]+=stmp; vs_unc[u] += "/";
											stmp.Form("%g",norm_up/norm_norminal);
											vs_unc[u]+=stmp; vs_unc[u] += "/";
											vs_unc[u]+=unc; vs_unc[u] += " ";
											if(debug>=10)cout<<vs_unc[u]<<endl;
										}
										//cout<<vs_unc[u]<<endl;
									}else{
										vs_unc[u]+=" - ";
									}
								}else{
									string tmp =GetUncertainy(c, t, vv_procnames, uncerlines[u]); 
									if(uncertypes[u]=="gmA" or uncertypes[u]=="gmN"){
										if(debug>=10)cout<<" binned histo: gamma "<<tmp<<endl;
										if(TString(tmp).IsFloat()){
											// make sure no negativae bin content
											float tmpd = TString(tmp).Atof()*(hnorm[t]->GetBinContent(r)<0?0:hnorm[t]->GetBinContent(r));
											TString tmps; 
											tmps.Form("%g",tmpd);
											tmp = tmps.Data();
										}
									}
									vs_unc[u]+=tmp; vs_unc[u] += " ";  // must be careful with gamma uncertainties
								}
								if(debug>=100) cout<<"debug "<<uncertypes[u]<< " "<<u<<endl;
							}
						}



						proc++;
					}
				}
				delete h;
				for(int t=0; t<n_proc; t++) {delete hn[t]; delete hnorm[t]; }
				/*
				   for(int t=0; t<n_proc; t++){
				   for(int u=0; u<nsyssources; u++){
				   if(hunc_up[t][u]) delete hunc_up[t][u];
				   if(hunc_up_norm[t][u]) delete hunc_up_norm[t][u];
				   if(hunc_dn[t][u]) delete hunc_dn[t][u];
				   if(hunc_dn_norm[t][u]) delete hunc_dn_norm[t][u];
				   }
				   }
				 */

			}// shape channels
			if(debug) cout<<" middle "<<endl;
			if(isShapeChannel==false) {
				n++;
				for(int t=0; t<vv_procnames[c].size(); t++){
					int totproc = GetTotProc(c, t, vv_procnames);
					s_bin+=channelnames[c]; s_bin+=" ";
					s_rate+=ss_old_rate[totproc+1];s_rate+=" ";
					// fill  s_process  //FIXME
					s_process += (t-nsigproc[c]+1) ; s_process += " ";
					s_process_name += vv_procnames[c][t] ; s_process_name += " ";

					for(int u=0; u<nsyssources; u++){
						if(uncertypes[u]=="param" || uncertypes[u]=="flatParam") continue;
						vector<string>	ss1 = uncerlines[u]; 
						if(ss1[0]!=uncernames[u]) continue;
						vs_unc[u]+=GetUncertainy(c, t, vv_procnames, ss1); vs_unc[u] += " ";  // must be careful with gamma uncertainties
					}
				}
			}
			if(debug) cout<<" done channel "<<channelnames[c]<<endl;
		}
		newlines.push_back(s_bin);
		newlines.push_back(s_process_name);
		newlines.push_back(s_process);
		newlines.push_back(s_rate);

		// extend uncertainty lines
		for(int u=0; u<nsyssources; u++){
			// line start with int number  "u"		
			newlines.push_back(vs_unc[u]);
		}

		for(int j=0; j<newlines.size(); j++){
			cardExpanded+=newlines[j]; 
			if(j<(newlines.size()-1)) cardExpanded+="\n";
		}
		if(debug){
			cout<<"***********print out of the new huge table****************"<<endl<<endl;
			cout<<cardExpanded<<endl<<endl<<endl;
		}
		ConfigureShapeModel(cms, mass, cardExpanded, histShapeLines, shapeuncertainties, parametricShapeLines, uncerlinesAffectingShapes, debug);
	}// line begin with "shape"
}

bool ConfigureModel(CountingModel *cms, double mass, const char* fileName, bool bUseHist, int debug){
      TString smass = ""; if(fmod(mass,1)==0) smass.Form("%.0f",mass);else smass.Form("%.1f",mass); 
	TString s = ReadFile(fileName);
	cout<<"data card name = "<<fileName<<endl;
	cms->SetModelName(fileName);
	return CheckIfDoingShapeAnalysis(cms, mass,  s, bUseHist, debug);
}


vector<TString> SplitIntoLines(TString ifileContentStripped, bool debug){	
	vector<TString> lines;
	//lines_array = ifileContentStripped.Tokenize(";");
	TObjArray *lines_array = ifileContentStripped.Tokenize("\n");
	TIterator* lineIt=lines_array->MakeIterator();
	TObject* line_o;

	const int nNeutrals=2;
	TString neutrals[nNeutrals]={"\t"," "};

	while((line_o=(*lineIt)())){

		TString line = (static_cast<TObjString*>(line_o))->GetString();

		// Strip spaces at the beginning and the end of the line
		line.Strip(TString::kBoth,' ');
		line.Strip(TString::kBoth,'\t');

		// Put the single statement in one single line
		line.ReplaceAll("\n","");
		line.ReplaceAll("\t"," ");

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

	return lines;
}
bool ConfigureShapeModel(CountingModel *cms, double mass, TString ifileContentStripped, 
		vector< vector<string> > histShapeLines, vector< vector<string> > histShapeUncLines,
		vector< vector<string> > parametricShapeLines, vector< vector<string> > uncerlinesAffectingShapes, int debug){
      TString smass = ""; if(fmod(mass,1)==0) smass.Form("%.0f",mass);else smass.Form("%.1f",mass); 
	// channel index starts from 1
	// systematics source index also starts from 1
	if(debug)for(int i=0; i<parametricShapeLines.size(); i++){
		for(int j=0; j<parametricShapeLines[i].size(); j++){
			cout<<parametricShapeLines[i][j]<<" ";
		}
		cout<<endl;
	}
	printf("\n *** \n");
	if(debug)for(int i=0; i<histShapeLines.size(); i++){
		for(int j=0; j<histShapeLines[i].size(); j++){
			cout<<histShapeLines[i][j]<<" ";
		}
		cout<<endl;
	}
	bool bUseHist = false;
	if(histShapeLines.size()){
		bUseHist=true;
		printf("\n *** \n");
		if(debug)for(int i=0; i<histShapeUncLines.size(); i++){
			for(int j=0; j<histShapeUncLines[i].size(); j++){
				cout<<histShapeUncLines[i][j]<<" ";
			}
			cout<<endl;
		}
	}


	TString duplicatingCard;
	vector<TString> duplicatingLines;

	//bool debug = true;
	// Now proceed with the parsing of the stripped file

	vector<TString> lines;
	lines = SplitIntoLines(ifileContentStripped, debug);
	//   
	//  For Andrey's input format 

	// get number of channels
	int nchannel = -1;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("imax ")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(debug)cout<<"NChannel = "<<ss[1]<<endl;
			TString s = ss[1];
			if(s.IsDigit()==false) {
				//cout<<"need a number after imax"<<", currently it's not a number but \'"<<s<<"\'"<<endl;
				//exit(0);
				nchannel = -1;
				cout<<"WARNING: You input of imax (number of channels) is not digit, hence LandS will not check the consistency"<<endl;
			}else {
				nchannel = s.Atoi();
			}

			duplicatingLines.push_back(lines[l]);
		} 
	}	
	if(nchannel==0) {
		cout<<"number of channels = 0,  please check imax"<<endl;
		return false;
	}

	// get observed dataset and cout how many channels in your model
	//
	vector<double> observeddata;
	bool hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("Observation ") or lines[l].BeginsWith("observation ")){
			observeddata.clear();
			if(hasFilled) cout<<"WARNING: You have two lines started with \"observation\", we will use the second line"<<endl;
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			for(int i=1; i<ss.size(); i++){
				if(TString(ss[i]).IsFloat()){
					double ev = (TString(ss[i])).Atof();
					observeddata.push_back(ev);
				}
			}
			hasFilled =  true;
			duplicatingLines.push_back(lines[l]);
		}
	}
	if(hasFilled==false) {cout<<"ERROR: need a line starting with \"Observation\" "<<endl; exit(0);}
	if(nchannel<0) nchannel = observeddata.size(); // if you don't assign a number to imax
	if(nchannel!= observeddata.size()) {
		cout<<"imax = "<<nchannel<<endl;
		cout<<"observeddata.size = "<<observeddata.size()<<endl;
		cout<<"ERROR: number of channels spcified in \"imax\" line is not consistent with \"observation\" line "<<endl;
		exit(0);
	}
	if(debug){
		cout<<"observed data: ";
		for(int i=0; i<nchannel; i++) cout<<observeddata[i]<<" ";
		cout<<endl;
	}

	// there must be one line of "process", which contains enumeration of processes in each channel
	// you might have another line with "process" which contains process names in each channel
	// since in one of the "process" lines,  for each channel, there must be one and only one process enumerated as "0"
	// we will count how many "0" in that line to determine how many channels actually go into your model
	// double check ....
	//
	// first check how many lines started with "process"
	int nlines_with_process = 0;
	vector< vector<string> > vss_processes; vss_processes.clear();
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("process ")){
			nlines_with_process++;
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			vss_processes.push_back(ss);
			duplicatingLines.push_back(lines[l]);
		}
	}
	if(nlines_with_process<1 or nlines_with_process>2 ){ cout<<"ERROR: you have "<<nlines_with_process<<" lines started with \"process\".  must be one or two"<<endl; exit(0); }
	// read number of channels from process line
	int tmpn = 0;  int l_proc = 0;
	vector<int> nsigproc; nsigproc.clear(); // record number of signal processes in each channel  (with process number <= 0)
	vector<string> processnames; processnames.clear();
	for(int l=0; l<vss_processes.size(); l++){
		for(int i=1; i<vss_processes[l].size(); i++){ // vss_processes[0] is "process"
			if(!TString(vss_processes[l][i]).IsFloat()) break; // you encounter comments or non-numbers, stop processing 
			if(tmpn==nsigproc.size())nsigproc.push_back(0);
			if(TString(vss_processes[l][i]).Atoi()<=0) nsigproc[nsigproc.size()-1]++;
			if(TString(vss_processes[l][i]).Atoi()==0) tmpn++; // as each channel has "0" standing for one of signal processes
			processnames.push_back(vss_processes[l][i]); // we temporarily assign process name as the "number", if another "process" line found,  we'll make change later
		}
		// debug
		for(int i=0; i<nsigproc.size(); i++){
			if(debug)cout<<" n signal process in channel "<<i<<": "<<nsigproc[i]<<endl;
		}


		if(tmpn>0) {
			l_proc = l; // record which line of "process" lines is for enumeration
			break; // already got the wanted line 
		}
	}
	if(nchannel != tmpn ) { cout<<"ERROR: \"observation\" and \"process\" are not consistent"<<endl; exit(0);}
	// read process names 
	if(vss_processes.size()==2){
		for(int i=1; i<vss_processes[1-l_proc].size(); i++){
			if(i>processnames.size()) break;
			processnames[i-1] = vss_processes[1-l_proc][i];
		}
	}

	// get numbers of processes in each channel
	int *nprocesses = new int[nchannel];
	for(int i=0; i<nchannel; i++) nprocesses[i]=0;
	bool hasLineStartWithBin = false;
	vector<string> channelnames; channelnames.clear();
	for(int l=0; l<lines.size(); l++){
		//if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==true){
		//	cout<<"WARNING:  there are two lines beginning with \"bin\" in your card. We will take the first one"<<endl;
		//}
		if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==false ){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size() < 1+2*nchannel) continue; // this "bin" line is too short, will process in below code

			// we now can proceed with names of bins instead of just numbers 1, 2, 3, ... N 
			int bin = 0;
			for(int i=1; i<ss.size(); i++){
				if(ss[i]!=ss[i-1]) { bin++; channelnames.push_back(ss[i]);	}
				if(bin<(nchannel+1) && bin>=1) nprocesses[bin-1]++;
			}
			if(bin != nchannel) {cout<<"ERROR: number of channels got from \"bin\" line is not consistent with \"observation\" line"<<endl; exit(0);};

			hasLineStartWithBin = true;
			duplicatingLines.push_back(lines[l]);
		} 
	}	
	if(!hasLineStartWithBin) {
		cout<<"Line beginning with \"bin\" is not found in your card"<<endl;
		exit(0);
	}
	if(debug){
		for(int i = 0; i<nchannel; i++) cout<<"processes in Channel "<<i<<" = "<<nprocesses[i]<<endl;
	}
	for(int i=0; i<nchannel; i++) if(nprocesses[i]==0) {cout<<"ERROR: channel "<<i<<", number of processes = 0"<<endl; exit(0);}

	// if you have "bins" or "binname" line,  then replace channel name  with the ones from this line
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("bins ") or lines[l].BeginsWith("binname") or lines[l].BeginsWith("bin ") ){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size() >= 1+nchannel*2) continue;// this "bin" line is too long, processed in above code 
			for(int i=1; i<ss.size(); i++){
				if(i>nchannel) break;
				channelnames[i-1]=ss[i];
			}
		}
	}

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
	hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("rate ")){
			TString tmps = TString::Format("%10s ","rate");
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size()-1<ntotprocesses) {cout<<"number of eventrate < "<<ntotprocesses<<endl; return false;} 
			for(int i=1; i<(ntotprocesses+1); i++){
				float ev = (TString(ss[i])).Atof();
				eventrate[i-1] = ev;
				tmps+= TString::Format("%7.2f ",ev);
			}
			hasFilled =  true;
			duplicatingLines.push_back(tmps);
		}
	}
	if(hasFilled==false) {cout<<"need a line starting with \"rate\" "<<endl; return false;}
	if(debug){
		cout<<"event rate: ";
		for(int i=0; i<ntotprocesses; i++) cout<<eventrate[i]<<" ";
		cout<<endl;
	}


	// FIXME we won't need this "kmax" keyword if we don't want to do a sanity check  ....
	//  get number of independant systematics sources
	int kmax = -1;
	//hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("kmax ")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size()>=2 && TString(ss[1]).IsFloat()){
				kmax = (TString(ss[1])).Atoi();
			}else{
				cout<<"WARNING: You input of kmax (number of systematic lines ) is not digit, hence LandS will not check the consistency"<<endl;
				kmax = -1;
			}
			//hasFilled = true;
			duplicatingLines.push_back(lines[l]);
		}
	}
	//if (hasFilled == false) cout<<"need \"kmax\" to set number of independant systematic sources "<<endl;
	if (kmax==0) cout<<"no systematics input, kmax=0"<<endl;

	if(debug) cout<<"number of independant systematics sources = "<<kmax<<endl;

	int INDEXofProcess = 1, INDEXofChannel=2;

	//  Fill the model with the yields
	index=0;
	for(int c=0; c<nchannel; c++){
		vector<double> sigbkgs; sigbkgs.clear();
		vector<string> tmpprocn; tmpprocn.clear();
		for(int p=0; p<nprocesses[c]; p++) {
			sigbkgs.push_back(eventrate[index]);
			tmpprocn.push_back(processnames[index]);
			index++;
		}
		if(debug) cout<<" nsigproc["<<c<<"] = "<<nsigproc[c]<<endl;

		bool isParametricChannel = false;
		bool isHistChannel = false;
		for(int i=0; i<parametricShapeLines.size(); i++){
			if(parametricShapeLines[i][2]==channelnames[c]){
				isParametricChannel = true;
			}
		}
		for(int i=0; i<histShapeLines.size(); i++){
			if(histShapeLines[i][2]==channelnames[c]){
				isHistChannel = true;
			}
		}
		if(isParametricChannel and isHistChannel) {cout<<"ERROR: "<<channelnames[c]<<" is both parametric and histogram channel, exit"<<endl; exit(1);}

		if(debug)cout<<" channel ["<<channelnames[c]<<"] isParametricChannel? "<<isParametricChannel<<endl;
		if(!isParametricChannel and !isHistChannel){
			cms->AddChannel(channelnames[c], sigbkgs, nsigproc[c]);
			cms->SetProcessNames(channelnames[c], tmpprocn);
			cms->AddObservedData(channelnames[c], observeddata[c]);
		}else if(isParametricChannel){
			if(TString(channelnames[c].substr(0,1)).IsAlpha()==false) {
				cout<<endl<<"ERROR: channel name ["<<channelnames[c]<<"] has a problem: "<<endl;
				cout<<"Parametric shape channel name should begin with alphabetai, exit"<<endl;
				exit(1);
			}
			vector<RooAbsPdf*> vspdf, vbpdf; vspdf.clear(); vbpdf.clear();
			vector<RooAbsArg*> vsExtraNorm, vbExtraNorm; vsExtraNorm.clear(); vbExtraNorm.clear();
			for(int i=0; i<sigbkgs.size(); i++){
				RooAbsPdf * pdf = (RooAbsPdf*)GetPdf(channelnames[c], tmpprocn[i], parametricShapeLines, mass);
				RooAbsArg * extraNom = (RooAbsArg*)GetExtraNorm(channelnames[c], tmpprocn[i], parametricShapeLines, mass);
				//pdf->SetName(TString(channelnames[c])+"_"+tmpprocn[i]+"_"+pdf->GetName()); // cause problem when 2l2q ggH and VBF use the same pdf 'signal'
				//if(extraNom)extraNom->SetName(TString(pdf->GetName())+"_"+"extranorm");
				//if(extraNom) extraNom->SetName(TString(channelnames[c])+extraNom->GetName());
				if(i<nsigproc[c]){
					vspdf.push_back(pdf);
					vsExtraNorm.push_back(extraNom);
				}else{
					vbpdf.push_back(pdf);
					vbExtraNorm.push_back(extraNom);
				}
			}
			vector<double> vsnorm, vbnorm; vsnorm.clear(); vbnorm.clear();
			for(int i=0; i<sigbkgs.size(); i++){
				if(i<nsigproc[c]) vsnorm.push_back(sigbkgs[i]);
				else vbnorm.push_back(sigbkgs[i]);
			}
			for(int i=0; i<vsnorm.size(); i++){
				if(debug)cout<<"* pdf "<<vspdf[i]->GetName()<<endl;
				if(debug)cout<<"* Normalization in data card = "<<vsnorm[i]<<endl;
				if(vsExtraNorm[i]) {
					//vsnorm[i]*=((RooAbsReal*)vsExtraNorm[i])->getVal();
					if(debug)cout<<"* ExtraNormalization in workspace = "<<((RooAbsReal*)vsExtraNorm[i])->getVal()<<endl;
				}
			}
			for(int i=0; i<vbnorm.size(); i++){
				if(debug)cout<<"* pdf "<<vbpdf[i]->GetName()<<endl;
				if(debug)cout<<"* Normalization in data card = "<<vbnorm[i]<<endl;
				if(vbExtraNorm[i]) {
				       	//vbnorm[i]*=((RooAbsReal*)vbExtraNorm[i])->getVal();
					if(debug)cout<<"* ExtraNormalization in workspace = "<<((RooAbsReal*)vbExtraNorm[i])->getVal()<<endl;
				}
			}
			RooAbsData *data = (RooAbsData*)GetRooAbsData(channelnames[c], "data_obs", parametricShapeLines, mass);
			if(data->sumEntries()!=observeddata[c] && observeddata[c]>=0) {
				cout<<"ERROR: In Channel: "<<channelnames[c]<<endl;
				cout<<"Observed data in card: "<<observeddata[c]<<" != "<<" RooAbsData.sumEntries"<<endl;
				exit(1);
			}

			cout<<" pdf observable ********* =  "<<endl;
			(vspdf[0]->getObservables(data))->Print();
			TString cuts;
			RooArgSet* observables = vspdf[0]->getObservables(data);
			std::auto_ptr<TIterator> iter(observables->createIterator());
			int ii = 0;
			for(RooRealVar *obs = (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
				cout<<" DELETEME *******  obs "<<obs->GetName()<<" bins = "<<obs->getBins()<<endl;
				//if(obs->getBins()>200) obs->setBins(200);
				if(ii>0) cuts+=" && ";
				cuts+=obs->GetName(); cuts+=">="; cuts+=obs->getMin(); cuts+=" && ";
				cuts+=obs->GetName(); cuts+="<="; cuts+=obs->getMax();
				ii++;
			}
			
			RooAbsData* reducedData = (RooAbsData*)data->reduce(*(vspdf[0]->getObservables(data)), cuts);
			//RooRealVar* x = dynamic_cast<RooRealVar*> (data->get()->first());
			RooRealVar* x = dynamic_cast<RooRealVar*> (observables->first());
			cms->AddChannel(channelnames[c], tmpprocn, x, vspdf, vsnorm, vsExtraNorm, vbpdf, vbnorm, vbExtraNorm);
			cms->AddObservedDataSet(channelnames[c], reducedData);
		}else if(isHistChannel){
			vector< vector<string> > shape_ch; //, uncer_ch;	
			if(TString(channelnames[c].substr(0,1)).IsAlpha()==false) {
				cout<<endl<<"ERROR: channel name ["<<channelnames[c]<<"] has a problem: "<<endl;
				cout<<"Histogram shape channel name should begin with alphabetai, exit"<<endl;
				exit(1);
			}

			for(int p=0; p<histShapeLines.size(); p++){
				if(histShapeLines[p][INDEXofChannel]!=channelnames[c] || histShapeLines[p][INDEXofProcess]!="data_obs") continue; // find the signal one
				if(TString(histShapeLines[p][4]).Contains(":")) continue; // parametric channel, will process later

				shape_ch.clear();
				for(int q=0; q<histShapeLines.size(); q++){
					if(histShapeLines[q][INDEXofChannel]==channelnames[c]) shape_ch.push_back(histShapeLines[q]);
				}
				//uncer_ch.clear();
				//for(int q=0; q<shapeuncertainties.size(); q++){
				//	if(shapeuncertainties[q][INDEXofChannel]==channelnames[c]) uncer_ch.push_back(shapeuncertainties[q]);
				//}
			}

			vector<TH1D*> vspdf, vbpdf; vspdf.clear(); vbpdf.clear();
			for(int i=0; i<sigbkgs.size(); i++){
				for(int q=0; q<shape_ch.size(); q++){
					if(shape_ch[q][INDEXofProcess]!=tmpprocn[i]) continue; 
					TH1D* pdf = (TH1D*)GetTObject(shape_ch[q][3], shape_ch[q][4]);

					if(i<nsigproc[c]){
						vspdf.push_back(pdf);
					}else{
						vbpdf.push_back(pdf);
					}
					if(pdf->Integral() == 0) { 
						cout<<"WARNING: channel ["<<channelnames[c]<<"] process ["<<tmpprocn[i]
							<<"] is shape, but the norminal histogram->Integral = 0"<<endl;
						//exit(0);
					}

					// check if normalization consistent with the rate declared in the data card
					double tmp_rate = sigbkgs[i];
					if(tmp_rate<0) {} //don't check normalization of histogram,  users' responsibility 
					else{ 
						if(fabs(tmp_rate-pdf->Integral())/pdf->Integral() > 1){ // if difference > 100%, then put Error and exit
							cout<<"ERROR: in channel ["<<channelnames[c]<<"] process ["<<tmpprocn[i]<<"], ";
							cout<<" declared rate in data card = "<< tmp_rate<<", not consistent with the integral of input histogram,";
							cout<<" which is "<<pdf->Integral()<<endl;
							cout<<" NOTE: if you don't want to check the consistency, you can put \"-\" in the datacard for the rate. "<<endl;
							cout<<"       In that case, you have the responsibility to make sure the integral (normalization) of histogram makes sense."<<endl;
							exit(0);
						}
					}
				}

			}	
			for(int q=0; q<shape_ch.size(); q++){
				if(shape_ch[q][INDEXofProcess]!="data_obs") continue; 
				TH1D* data= (TH1D*)GetTObject(shape_ch[q][3], shape_ch[q][4]);
				if(observeddata[c]>=0 and data->Integral()!=observeddata[c]) {
					cout<<"ERROR: "<<channelnames[c]<<" data in cards "<<observeddata[c]<<" != hist integral "<<data->Integral()<<endl;
					cout<<" you can set data to be negative number in the text card, so that the tool won't check and will take the data hist directly. exit"<<endl;
					exit(1);
				}
				cms->AddChannel(channelnames[c], tmpprocn, vspdf, vbpdf);
				cms->SetProcessNamesTH(channelnames[c], tmpprocn);
				cms->AddObservedDataTH(channelnames[c], data);
				break;
			}
		}
	}	

	//cms->Print();

	if(debug) cout<<"filled yields"<<endl;
	// fill the model with systemics

	// sections are separated by "---------"  in the data card
	// last section are all uncertainties,  it's analyzers' responsibility to make it right  
	int nsyssources = 0;
	for(int j=lines.size()-1; j>=0; j--){
		if(lines[j].BeginsWith("---") )	break;
		if(lines[j].BeginsWith("rate "))	break;
		if(lines[j].BeginsWith("bin ") )	break;
		if(lines[j].BeginsWith("process "))break;
		if(lines[j].BeginsWith("imax "))break;
		if(lines[j].BeginsWith("jmax "))break;
		if(lines[j].BeginsWith("kmax "))break;
		nsyssources++;
	}
	if(kmax>=0 && kmax!=nsyssources) {cout<<"kmax !=  number of independant uncertainties"<<endl; exit(0);}

	vector< vector<string> > vshape_params_unclines; vshape_params_unclines.clear();
	if(debug)cout<<" total "<<nsyssources<<" systs"<<endl;
	for(int s=0; s<nsyssources; s++){
		vector<string>	ss; 

		int hasUncSourceWithIndex_s = -1;
		hasUncSourceWithIndex_s = lines.size() - nsyssources + s;
		ss.clear();
		StringSplit(ss, lines[hasUncSourceWithIndex_s].Data(), " ");
		if(debug)cout<<" @@@@@@@ systemticline: "<<lines[hasUncSourceWithIndex_s]<<endl;

		TString tmps; tmps.Form("%s ", ss[0].c_str());

		string indexcorrelation = ss[0];
		TString sTmp=indexcorrelation;
		bool bUncIsFloatInFit = true;
		if(sTmp.Contains("[nofloat]")) bUncIsFloatInFit=false;
		indexcorrelation = sTmp.ReplaceAll("[nofloat]","").Data();

		bool bHistShapeUnc = false;

		int pdf = 0; 
		// see CountingModel.h
		if(ss[1]=="lnN") pdf=typeLogNormal;   // typeLogNormal
		else if(ss[1]=="flat") pdf=typeFlat; // typeFlat
		else if(ss[1]=="trG") pdf=typeTruncatedGaussian; // typeTruncatedGaussian
		else if(ss[1]=="gmA" or ss[1]=="gmN" or ss[1]=="gmM") pdf=typeGamma; // typeControlSampleInferredLogNormal;  gmA was chosen randomly, while in gmN, N stands for yield in control sample;  gmM stands for Multiplicative gamma distribution, it implies using a Gamma distribution not for a yield but for a multiplicative correction
		else if(ss[1]=="shapeN" or ss[1]=="shapeN2" or ss[1]=="shapeStat") { pdf=typeLogNormal; bHistShapeUnc = (bUseHist?true:false); } 
		else if(ss[1]=="shape" or ss[1]=="shape2" or ss[1]=="shapeQ" or ss[1]=="shape?") { pdf=typeShapeGaussianQuadraticMorph; bHistShapeUnc=(bUseHist?true:false);}
		else if(ss[1]=="shapeL" ) { pdf=typeShapeGaussianLinearMorph; bHistShapeUnc= (bUseHist?true:false);}
		else if(ss[1]=="param" ) {
			pdf=typeBifurcatedGaussian;
			vshape_params_unclines.push_back(ss);
		}else if(ss[1]=="flatParam") {
			pdf=typeFlat;
			vshape_params_unclines.push_back(ss);
		}
		else {
			//pdf =  (TString(ss[1])).Atoi();
			cout<<" ERROR: Uncertainty type ["<<ss[1]<<"] is not defined yet "<<endl;
			exit(1);
		}

		if(debug) cout<<" pdf type "<<pdf<<endl;

		tmps+= TString::Format("%8s ", ss[1].c_str());
		if(pdf==typeGamma && ss[1]!="gmM") {
			tmps+= TString::Format("%8s ", ss[2].c_str());
			if(ss.size()<ntotprocesses+3) { cout<<"Error... uncertainty "<<ss[0]<<": doesn't have enough collums"<<endl; exit(1) ; }
		}
		else if(pdf==typeBifurcatedGaussian || pdf==typeFlat) {
			for(int ii=2; ii<ss.size(); ii++) {tmps+=ss[2]; tmps+=" ";}
		}
		else		   tmps+= "         ";

		if(ss[1]=="flat" or ss[1]=="lnN" or ss[1]=="trG" or ss[1]=="gmM" or ss[1]=="shapeN" or ss[1]=="shapeN2" or ss[1]=="shape" or ss[1]=="shape2" or ss[1]=="shapeStat" or ss[1]=="shapeL" 
				or ss[1]=="shapeQ" or ss[1]=="shape?"){
			if(ss.size()<ntotprocesses+2) { cout<<"Error... uncertainty "<<ss[0]<<": doesn't have enough collums"<<endl; exit(1) ; }
		}

		bool filledThisSource = false;

		int pdf_old = pdf;
		for(int p=0; p<ntotprocesses; p++){
			pdf=pdf_old;

			double err, errup, rho;  // to allow asymetric uncertainties
			double shape[8];

			string channelName = channelnames[binnumber[p]-1];
			bool isParametricChannel = false;
			bool isHistChannel = false;
			for(int i=0; i<parametricShapeLines.size(); i++){
				if(parametricShapeLines[i][2]==channelName){
					isParametricChannel = true;
				}
			}
			for(int i=0; i<histShapeLines.size(); i++){
				if(histShapeLines[i][2]==channelName){
					isHistChannel = true;
				}
			}

			if(pdf==typeLogNormal){
				if(ss[p+2]=="-" || TString(ss[p+2]).Atof()==0) {
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/")){
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, ss[p+2], "/");
					if((TString(asymetricerrors[0])).Atof()<=0) { //cout<<"ERROR:  Kappa can't be <=0 "<<endl; exit(0); ;
						err = 1;
					}else	err= 1./(TString(asymetricerrors[0])).Atof()-1.0; // downside 
					
					if((TString(asymetricerrors[0])).Atof()<=0) { //cout<<"ERROR:  Kappa can't be <=0 "<<endl; exit(0); ;
						errup =1;
					}else errup= (TString(asymetricerrors[1])).Atof()-1.0;  // upside
				}else {
					err= (TString(ss[p+2])).Atof()-1.0;
					errup = err;
				}
				tmps+= TString::Format("%7.2f ",err+1.0);
				if(err <= -1 or errup <= -1) {
					cout<<"ERROR: Kappa in lognormal  can't be negative, please check your input card at "<<s+1<<"th source, "<<p+2<<"th entry"<<endl;
					exit(0);
				}
				if(err == 0. && errup==0.) continue;
			}
			else if(pdf==typeFlat && ss[1]=="flat"){
				if(ss[p+2]=="-") {
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/")){
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, ss[p+2], "/");
					if((TString(asymetricerrors[0])).Atof()<0) {
						err = 1;
					}else
						err= (TString(asymetricerrors[0])).Atof(); // downside 
					if((TString(asymetricerrors[1])).Atof()<0) {
						errup = 1;
					}else
						errup= (TString(asymetricerrors[1])).Atof();  // upside
				}else {
					tmps+= TString(ss[p+2]);
					cout<<"WARNING: flat uncertainty on normalization need to use \"/\" to separate min and max factors.  being skiped"<<endl;
					continue;
				}
				tmps+= TString(ss[p+2]);
			}
			else if(pdf==typeTruncatedGaussian){
				if(ss[p+2]=="-"){
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/")){
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, ss[p+2], "/");
					err= (TString(asymetricerrors[0])).Atof();
					errup= (TString(asymetricerrors[1])).Atof();
				}else {
					err= (TString(ss[p+2])).Atof();
					errup = err;
				}
				tmps+= TString::Format("%7.2f ",err);
				if(err == 0.) continue;
			}else if(pdf==typeShapeGaussianQuadraticMorph or pdf==typeShapeGaussianLinearMorph){
				if(ss[p+2]=="-"){
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/") ){
					// uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm, scale_down_of_gaussian,  siglebin_or_binned
					vector<string> params; params.clear();
					StringSplit(params, ss[p+2], "/");
					if(params.size()==7) {
						for(int i=0; i<7; i++) shape[i]=TString(params[i]).Atof();
						shape[7]=1.; // binned
					}else if(params.size()==2){
						pdf = typeLogNormal; // FIXME  if "shape" also affect single cut and count or parametric channel, then tranfer it to LogNormal 
						for(int i=0; i<6; i++) {
							if(i<4)shape[i]= 0 ;
							else shape[i]=TString(params[i-4]).Atof();
						}
						if((TString(params[0])).Atof()<=0) {//{ cout<<"ERROR:  Kappa can't be <=0 "<<":  "<<asymetricerrors[0]<<endl; exit(0); };
							err = 1;
						}else
							err = 1./TString(params[0]).Atof() -1.;
						if((TString(params[1])).Atof()<=0) {//{ cout<<"ERROR:  Kappa can't be <=0 "<<":  "<<asymetricerrors[0]<<endl; exit(0); };
							errup = 1;
						}else
							errup = TString(params[1]).Atof() -1.;

						shape[6]=0;
						shape[7]=0;// cut_and_count
					}else { cout<<"ERROR: typeShape input format not correct "<<endl; exit(0);}
				}else{
					for(int i=0; i<4; i++) {
						shape[i]= 0 ;
					}
					shape[4]=TString(ss[p+2]).Atof();
					shape[5]=shape[4];
					shape[6]=0;
					shape[7]=0;
					if(isHistChannel and bHistShapeUnc){}else pdf = typeLogNormal;
					err = TString(ss[p+2]).Atof()-1;
					errup = err;
				}
				tmps= TString::Format("%7.2f ", shape[0]);
			}
			else if(pdf==typeGamma){
				if(ss[1]=="gmA" or ss[1]=="gmN"){
					if(ss[p+3]=="-"){
						tmps+= "    -   ";
						continue;
					}
					rho= (TString(ss[p+3])).Atof(); // the factor between control region and signal region
					tmps+= TString::Format("%7.2f ",rho);
					if(rho < 0) {
						cout<<"Alpha in gamma can't be negative, please check your input card at "<<s+1<<"th source, "<<p+3<<"th entry"<<endl;
						exit(0);
					}
					if(rho == 0.) continue;
				}
				else if(ss[1]=="gmM"){
					if(ss[p+2]=="-") {
						tmps+= "    -   ";
						continue;
					}
					err= (TString(ss[p+2])).Atof();
					tmps+= TString::Format("%7.2f ",err);
					if(err < 0 || err>1 ) {
						cout<<"relative error in gamma function can't be negative or > 100%, please check your input card at "<<s+1<<"th source, "<<p+2<<"th entry"<<endl;
						exit(0);
					}
					if(err == 0.) continue;
				}else{
					cout<<"Gamma function assumed,  the acronym should be gmA, gmN or gmM,  not "<< ss[1]<<", exit "<<endl;
					exit(0);
				}
			}

			//if(debug)cout<<" To add unc"<<endl;

			//cout<<" DELETEME 1 pdftype: "<<pdf<<endl;
			if( pdf==typeLogNormal ||pdf==typeTruncatedGaussian || (pdf==typeFlat && ss[1]=="flat")){
				if(isParametricChannel)cms->AddUncertaintyOnShapeNorm(channelName, subprocess[p], err, errup, pdf, indexcorrelation );
				else { 
					if(!isHistChannel){
						if(ss[1]=="shapeStat"){
							TString stmpunc = indexcorrelation; stmpunc+=channelName; stmpunc+=subprocess[p];
							cms->AddUncertainty(channelName, subprocess[p], err, errup, pdf, stmpunc.Data() );
							cms->TagUncertaintyFloatInFit(stmpunc.Data(), bUncIsFloatInFit);
						}else{	cms->AddUncertainty(channelName, subprocess[p], err, errup, pdf, indexcorrelation );
							cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
						}
					}else{
						TString sunc = ss[p+2];
						if(bHistShapeUnc and sunc.IsFloat() && sunc.Atof()>0){
							double scaleUnc  = sunc.Atof();
							bool filled = false;
							for(int i=0; i<histShapeUncLines.size(); i++){ // look for the histograms from the list 
								if( histShapeUncLines[i][INDEXofChannel]==channelName and
										histShapeUncLines[i][INDEXofProcess]==processnames[p]
										and histShapeUncLines[i][4]==ss[0]){
									TH1D * hunc_up= (TH1D*)GetTObjectA(histShapeUncLines[i][3], histShapeUncLines[i][5]);
									if(hunc_up==NULL) continue;

									TH1D *hn= cms->GetExpectedTH(channelName, processnames[p]);	
									if(hn==NULL) {
										cout<<" channle ["<<channelName<<"] proc ["<<processnames[p]<<"]"<<" histogram is not found "<<endl;
										exit(1);
									}
									TString sname  = "del_clone"; sname+=i;
									TH1D* hnorm=(TH1D*)hn->Clone(sname);
									if(hn->Integral()!=0)hnorm->Scale(1./hn->Integral());

									if(hunc_up->Integral()== 0 && hnorm->Integral()!=0) { 
										cout<<" hnorm "<<hnorm->GetName()<<": integral = "<<hnorm->Integral()<<endl;
										cout<<"ERROR: channel ["<<channelName<<"] process ["<<processnames[p]
											<<"] is shape, but the up_shift histogram->Integral = 0"<<endl;
										exit(1);
									}
									TH1D* hunc_dn= (TH1D*)GetTObjectA(histShapeUncLines[i][3], histShapeUncLines[i][6]);
									if(hunc_dn==NULL) continue;
									if(hunc_dn->Integral()== 0 && hnorm->Integral()!=0) { 
										cout<<"ERROR: channel ["<<channelName<<"] process ["<<processnames[p]
											<<"] is shape, but the down_shift histogram->Integral = 0"<<endl;
										exit(1);
									}
									sname  = "hup_clone"; sname+=i; sname+=p;
									TH1D* hunc_up_norm=(TH1D*)hunc_up->Clone(sname);
									if(hunc_up->Integral()!=0)hunc_up_norm->Scale(1./hunc_up->Integral());
									sname  = "hdn_clone"; sname+=i; sname+=p;
									TH1D* hunc_dn_norm=(TH1D*)hunc_dn->Clone(sname);
									if(hunc_dn->Integral()!=0)hunc_dn_norm->Scale(1./hunc_dn->Integral());

									if(hnorm->Integral()==0){
										cout<<" DEBUGGGGGGG 2"<<endl;
										for(int ii=0; ii<=hunc_up_norm->GetNbinsX(); ii++){
											hunc_dn_norm->SetBinContent(ii,0);
											hunc_up_norm->SetBinContent(ii,0);
										}
									}
									if(0&&debug){
										cout<<"channel ["<<channelName<<"] process ["<<processnames[p]<<": sys."<<ss[0]<<endl;
										cout<<hnorm->GetName()<<":"<<hnorm->Integral()<<endl;
										cout<<hunc_dn_norm->GetName()<<":"<<hunc_dn_norm->Integral()<<endl;
										cout<<hunc_up_norm->GetName()<<":"<<hunc_up_norm->Integral()<<endl;
									}
									if(filled) {cout<<" filling twice c:"<<channelName<<" p:"<<processnames[p]<<" sys:"<<ss[0]<<endl; exit(1);}
									filled = true;

									// shapeN, shapeN2, or shapeStat(?) 
									for(int b=1;b<=hunc_up_norm->GetNbinsX(); b++){
										hunc_up_norm->SetBinContent(b, hnorm->GetBinContent(b)==0?1:hunc_up_norm->GetBinContent(b)/hnorm->GetBinContent(b));
										hunc_dn_norm->SetBinContent(b, hnorm->GetBinContent(b)==0?1:hunc_dn_norm->GetBinContent(b)/hnorm->GetBinContent(b));
									}


									vector<TH1D*> vunc; vunc.clear();
									vunc.push_back(hunc_dn_norm);
									vunc.push_back(hunc_up_norm);

									//vunc.push_back(hnorm);
									//TH1D* htmp = new TH1D("htmp","htmp", 1,0,1);
									//htmp->SetBinContent(1,hn->Integral());
									//vunc.push_back((TH1D*)htmp->Clone());
									//htmp->SetBinContent(1,hunc_dn->Integral()/hn->Integral());
									//vunc.push_back((TH1D*)htmp->Clone());
									//htmp->SetBinContent(1,hunc_up->Integral()/hn->Integral());
									//vunc.push_back((TH1D*)htmp->Clone());
									//htmp->SetBinContent(1,scaleUnc);
									//vunc.push_back((TH1D*)htmp->Clone());

									cms->AddUncertaintyTH(channelName, subprocess[p], vunc, pdf, indexcorrelation);
									cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
								}
							}	
							if(!filled and !TString(ss[1]).Contains("?")) {
								cout<<"ERROR channel ["<<channelName<<"] process ["<<processnames[p]<<"] sys ["<<ss[0]
									<<"] is shape uncertainty, we need corresponding histogram "<<endl;
								exit(1);
							}
						}else{
							cms->AddUncertaintyTH(channelName, subprocess[p], err, errup, pdf, indexcorrelation );
							cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
						}
					}
				}
				filledThisSource = true;
			}
			if(pdf==typeGamma){
				if(ss[1]=="gmA" or ss[1]=="gmN"){
					//if(debug)cout<<"  it's  gmN "<<endl;
					double N = (TString(ss[2]).Atof()); // your input should be Most Probable Value,  while mean value is N+1
					if(N<0)  {
						cout<<"Yield in control Region in gamma can't be negative, please check your input card at "<<s+1<<"th source, "<<2<<"th entry"<<endl;
						exit(0);
					}
					//if(debug)cout<<" ch "<<channelName<<", subp  "<<p<<" "<<" rho="<<rho<<" N="<<N<<" indexcorrelation="<<indexcorrelation<<endl;
					if(isParametricChannel)cms->AddUncertaintyOnShapeNorm(channelName, subprocess[p], rho, 0, N, pdf, indexcorrelation );
					else if(isHistChannel){
						cms->AddUncertaintyTH(channelName, subprocess[p], rho, 0, N, pdf, indexcorrelation );
					}else 
						cms->AddUncertainty(channelName, subprocess[p], rho, 0, N, pdf, indexcorrelation );
					//FIXME  here need to check:  rho*N == cms->GetExpectedNumber(channelName, subprocess[p]);
				}
				else if(ss[1]=="gmM"){
					double N = 1./err/err-1; // we use convention: mean of events is 1./err/err,  so we do "-1" here, and then add back "1"  in src/CountingModel.cc
					//double N = (int) (1./err/err);
					if(isParametricChannel)cms->AddUncertaintyOnShapeNorm(channelName, subprocess[p], -1, 0, N, pdf, indexcorrelation ); // use "rho = -1" here to imply that this gamma term is not for control sample inferred uncertainties, but multiplicative gamma function ....
					else if(isHistChannel){
						cms->AddUncertaintyTH(channelName, subprocess[p], -1, 0, N, pdf, indexcorrelation ); 
					}else cms->AddUncertainty(channelName, subprocess[p], -1, 0, N, pdf, indexcorrelation ); // use "rho = -1" here to imply that this gamma term is not for control sample inferred uncertainties, but multiplicative gamma function ....
					//FIXME  here need to check:  all uncertainties with same indexcorrelation are the same,   can't be different currently ... 
					//  e.g.   check   N == cms->Get_v_Gamma(indexcorrelation); // can't do here before ConfigUncertaintyPdfs()
					//we might allow them different and do rescaling 
				}
				//if(debug)cout<<"  added  gamma "<<endl;
				filledThisSource = true;
			}
			if(pdf==typeShapeGaussianLinearMorph or pdf==typeShapeGaussianQuadraticMorph ){
				if(isParametricChannel) { cout<< "WARNING  Morph is not working for parametric shape input ,  will skip  this uncertainty" << endl; continue; }
				else{
					if(!bUseHist)cms->AddUncertainty(channelName, subprocess[p], 8, shape, pdf, indexcorrelation );
					else{
						TString sunc = ss[p+2];
						if(bHistShapeUnc and isHistChannel and sunc.IsFloat() && sunc.Atof()>0){
							double scaleUnc  = sunc.Atof();
							bool filled = false;
							for(int i=0; i<histShapeUncLines.size(); i++){ // look for the histograms from the list 
								if( histShapeUncLines[i][INDEXofChannel]==channelName and
										histShapeUncLines[i][INDEXofProcess]==processnames[p]
										and histShapeUncLines[i][4]==ss[0]){
									TH1D * hunc_up= (TH1D*)GetTObjectA(histShapeUncLines[i][3], histShapeUncLines[i][5]);
									if(hunc_up==NULL) continue;

									TH1D *hn= cms->GetExpectedTH(channelName, processnames[p]);	
									if(hn==NULL) {
										cout<<" channle ["<<channelName<<"] proc ["<<processnames[p]<<"]"<<" histogram is not found "<<endl;
										exit(1);
									}
									TString sname  = "ddel_clone"; sname+=i;
									TH1D* hnorm=(TH1D*)hn->Clone(sname);
									if(hn->Integral()!=0)hnorm->Scale(1./hn->Integral());

									if(hunc_up->Integral()== 0 && hnorm->Integral()!=0) { 
										cout<<" hnorm "<<hnorm->GetName()<<": integral = "<<hnorm->Integral()<<endl;
										cout<<"ERROR: channel ["<<channelName<<"] process ["<<processnames[p]
											<<"] is shape, but the up_shift histogram->Integral = 0"<<endl;
										exit(1);
									}
									TH1D* hunc_dn= (TH1D*)GetTObjectA(histShapeUncLines[i][3], histShapeUncLines[i][6]);
									if(hunc_dn==NULL) continue;
									if(hunc_dn->Integral()== 0 && hnorm->Integral()!=0) { 
										cout<<"ERROR: channel ["<<channelName<<"] process ["<<processnames[p]
											<<"] is shape, but the down_shift histogram->Integral = 0"<<endl;
										exit(1);
									}
									sname  = "hhup_clone"; sname+=i; sname+=p;
									TH1D* hunc_up_norm=(TH1D*)hunc_up->Clone(sname);
									if(hunc_up->Integral()!=0)hunc_up_norm->Scale(1./hunc_up->Integral());
									sname  = "hhdn_clone"; sname+=i; sname+=p;
									TH1D* hunc_dn_norm=(TH1D*)hunc_dn->Clone(sname);
									if(hunc_dn->Integral()!=0)hunc_dn_norm->Scale(1./hunc_dn->Integral());

									if(hnorm->Integral()==0){
										cout<<" DEBUGGGGGGG 2"<<endl;
										for(int ii=0; ii<=hunc_up_norm->GetNbinsX(); ii++){
											hunc_dn_norm->SetBinContent(ii,0);
											hunc_up_norm->SetBinContent(ii,0);
										}
									}
									if(0&&debug){
										cout<<"channel ["<<channelName<<"] process ["<<processnames[p]<<": sys."<<ss[0]<<endl;
										cout<<hnorm->GetName()<<":"<<hnorm->Integral()<<endl;
										cout<<hunc_dn_norm->GetName()<<":"<<hunc_dn_norm->Integral()<<endl;
										cout<<hunc_up_norm->GetName()<<":"<<hunc_up_norm->Integral()<<endl;
									}
									if(filled) {cout<<" filling twice c:"<<channelName<<" p:"<<processnames[p]<<" sys:"<<ss[0]<<endl; exit(1);}
									filled = true;

									// shapeN, shapeN2, or shapeStat(?) 
									//for(int b=1;b<=hunc_up_norm->GetNbinsX(); b++){
									//	hunc_up_norm->SetBinContent(b, hnorm->GetBinContent(b)==0?1:hunc_up_norm->GetBinContent(b)/hnorm->GetBinContent(b));
									//	hunc_dn_norm->SetBinContent(b, hnorm->GetBinContent(b)==0?1:hunc_dn_norm->GetBinContent(b)/hnorm->GetBinContent(b));
									//}
									
									//vs_unc[u]+=down; vs_unc[u] += "/";
									//stmp.Form("%g",up);
									//vs_unc[u]+=stmp; vs_unc[u] += "/";
									//stmp.Form("%g",norminal);
									//vs_unc[u]+=stmp; vs_unc[u] += "/";
									//stmp.Form("%g",norm_norminal);
									//vs_unc[u]+=stmp; vs_unc[u] += "/";
									//stmp.Form("%g",norm_down/norm_norminal);
									//vs_unc[u]+=stmp; vs_unc[u] += "/";
									//stmp.Form("%g",norm_up/norm_norminal);
									//vs_unc[u]+=stmp; vs_unc[u] += "/";
									//vs_unc[u]+=unc; vs_unc[u] += " ";
									
									// shape, shapeL, shapeQ
									vector<TH1D*> vunc; vunc.clear();
									vunc.push_back(hunc_dn_norm);
									vunc.push_back(hunc_up_norm);
									vunc.push_back(hnorm);
									TString stmp = channelName; stmp+="_"; stmp+=p;stmp+=s; stmp+=i;
									TH1D* htmp = new TH1D(stmp,"htmp", 1,0,1);
									htmp->SetBinContent(1,hn->Integral());
									stmp = "htmp"; stmp+=p; stmp+=s; stmp+=i; stmp+="_0";
									vunc.push_back((TH1D*)htmp->Clone(stmp));
									htmp->SetBinContent(1,hunc_dn->Integral()/hn->Integral());
									stmp = "htmp"; stmp+=p; stmp+=s; stmp+=i; stmp+="_1";
									vunc.push_back((TH1D*)htmp->Clone());
									htmp->SetBinContent(1,hunc_up->Integral()/hn->Integral());
									stmp = "htmp"; stmp+=p; stmp+=s; stmp+=i; stmp+="_2";
									vunc.push_back((TH1D*)htmp->Clone());
									htmp->SetBinContent(1,scaleUnc);
									stmp = "htmp"; stmp+=p; stmp+=s; stmp+=i; stmp+="_3";
									vunc.push_back((TH1D*)htmp->Clone());
									delete htmp;

									cms->AddUncertaintyTH(channelName, subprocess[p], vunc, pdf, indexcorrelation);
								}
							}	
							if(!filled and !TString(ss[1]).Contains("?")) {
								cout<<"ERROR channel ["<<channelName<<"] process ["<<processnames[p]<<"] sys ["<<ss[0]
									<<"] is shape uncertainty, we need corresponding histogram "<<endl;
								exit(1);
							}

						}
					}
				}
				filledThisSource = true;
			}

		}
		if(!filledThisSource) {
			if(debug)cout<<"WARNING: The "<< s+1 <<"th source of uncertainties are all 0. "<<endl;
		}else{
			if(pdf!=typeLogNormal)cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
		}

		if(pdf==typeBifurcatedGaussian){
			if(ss.size()<4) {
				cout<<"ERROR uncertainty line on parameter: "<<ss[0]<<" don't have enough parameters, need at least mean and sigma"<<endl;
				for(int ii=0; ii<ss.size(); ii++){
					cout<<ss[ii]<<" ";
				}
				cout<<endl;
				exit(0);
			}
			double mean = TString(ss[2]).Atof();
			double sigmaL, sigmaR;
			double rangeMin=0, rangeMax=0;
			if(TString(ss[3]).Contains("/")){
				vector<string> asymetricerrors; asymetricerrors.clear();
				StringSplit(asymetricerrors, ss[3], "/");
				sigmaL = (TString(asymetricerrors[0])).Atof(); // downside 
				sigmaR = (TString(asymetricerrors[1])).Atof();  // upside
			}else {
				sigmaL = (TString(ss[3])).Atof();
				sigmaR = sigmaL;
			}
			if(ss.size()>4){
				TString sr = ss[4];
				if(sr.BeginsWith("[") and sr.EndsWith("]") and sr.Contains(",") ){
					sr.ReplaceAll("[", "");
					sr.ReplaceAll("]", "");
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, sr.Data(), ",");
					rangeMin= (TString(asymetricerrors[0])).Atof(); // downside 
					rangeMax= (TString(asymetricerrors[1])).Atof();  // upside
				}
			}
			bool succes = cms->AddUncertaintyOnShapeParam(indexcorrelation, mean, sigmaL, sigmaR, rangeMin, rangeMax );
			if(succes)cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
		}
		if(pdf==typeFlat && ss[1]=="flatParam"){
			if(debug)cout<<" typeFlat "<<endl;
			// FIXME add hoc implementation of flat param  uncertainties ....  need to think about also for normalization terms
			bool succes= false;
			if(ss.size()<2) {cout<<"ERROR: parameter line with type flatParam should have at least two columns "<<endl; exit(1);}
			else if(ss.size()==2){
				succes=cms->AddUncertaintyOnShapeParam(indexcorrelation);
			}else if(ss.size()==3){
				double norminalValue = TString(ss[2]).Atof();
				succes=cms->AddFlatParam(indexcorrelation, norminalValue, norminalValue, norminalValue);// fix the param at norminalValue
			}else if(ss.size()>3){
				double norminalValue = TString(ss[2]).Atof();
				double rangeMin=norminalValue, rangeMax=norminalValue;
				TString sr = ss[3];
				if(sr.BeginsWith("[") and sr.EndsWith("]") and sr.Contains(",") ){
					sr.ReplaceAll("[", "");
					sr.ReplaceAll("]", "");
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, sr.Data(), ",");
					rangeMin= (TString(asymetricerrors[0])).Atof(); // downside 
					rangeMax= (TString(asymetricerrors[1])).Atof();  // upside
				}else{
					cout<<" Error:  flatParam format should be:  "<<endl;
					cout<<" UncName flatParam norminalValue [min,max] "<<endl;
					exit(1); 
				}
				succes=cms->AddFlatParam(indexcorrelation, norminalValue, rangeMin, rangeMax);// fix the param at norminalValue
			}
			if(succes)cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
		}

		duplicatingLines.push_back(tmps);
	}
	for(int i=0; i<uncerlinesAffectingShapes.size(); i++){
		vector<string>	ss = uncerlinesAffectingShapes[i]; 
		string indexcorrelation = ss[0];
		TString sTmp=indexcorrelation;
		bool bUncIsFloatInFit = true;
		if(sTmp.Contains("[nofloat]")) bUncIsFloatInFit=false;
		indexcorrelation = sTmp.ReplaceAll("[nofloat]","").Data();

		if(ss.size()<6) {
			cout<<"ERROR format of uncertainty affecting parameters should be as follows"<<endl;
			cout<<"<systname> affects <bin> <process> <parameter> <effect strength>"<<endl;
			cout<<"      while in data card, following line not compatible with the convention: "<<endl;
			for(int ii=0; ii<ss.size(); ii++){
				cout<<ss[ii]<<" ";
			}
			cout<<endl;
			exit(0);
		}
		//double mean = TString(ss[3]).Atof();
		double sigmaL, sigmaR;
		double rangeMin=0, rangeMax=0;
		if(TString(ss[5]).Contains("/")){
			vector<string> asymetricerrors; asymetricerrors.clear();
			StringSplit(asymetricerrors, ss[5], "/");
			sigmaL = (TString(asymetricerrors[0])).Atof(); // downside 
			sigmaR = (TString(asymetricerrors[1])).Atof();  // upside
		}else {
			sigmaL = (TString(ss[5])).Atof();
			sigmaR = sigmaL;
		}
		if(ss.size()>6){
			TString sr = ss[6];
			if(sr.BeginsWith("[") and sr.EndsWith("]") and sr.Contains(",") ){
				sr.ReplaceAll("[", "");
				sr.ReplaceAll("]", "");
				vector<string> asymetricerrors; asymetricerrors.clear();
				StringSplit(asymetricerrors, sr.Data(), ",");
				rangeMin= (TString(asymetricerrors[0])).Atof(); // downside 
				rangeMax= (TString(asymetricerrors[1])).Atof();  // upside
			}
		}
		//cms->AddUncertaintyAffectingShapeParam(ss[0], ss[2], mean, sigmaL, sigmaR, rangeMin, rangeMax );
		cms->AddUncertaintyAffectingShapeParam(indexcorrelation, ss[4], sigmaL, sigmaR);
		cms->TagUncertaintyFloatInFit(indexcorrelation, bUncIsFloatInFit);
	}

	if(debug) cout<<"filled systematics"<<endl;

	//cms->Print(100);
	if(debug) {
		cout<<"start duplicating this card:"<<endl<<endl;
		for(int i=0; i<duplicatingLines.size(); i++) cout<<duplicatingLines[i]<<endl;
		cout<<endl<<endl<<"end duplicating this card:"<<endl<<endl;
	}

	if(nprocesses) delete [] nprocesses;
	if(binnumber) delete [] binnumber;
	if(subprocess) delete [] subprocess;
	if(eventrate) delete [] eventrate;
	return true;
}
TString GetWordFromLine(TString line, int index, string delim){
	vector<string>	ss; 
	ss.clear();
	StringSplit(ss, line.Data(), delim);
	if(index>=ss.size()) {
		cout<<"ERROR You are trying to access the "<<index+1<<"th word from line below, which has only "<<ss.size()<<" words: "<<endl;
		cout<<line<<endl;
		exit(0);
	}
	TString s = ss[index];
	return s;
}
RooAbsPdf* GetPdf(string c, string p, const vector< vector<string> >& lines, double mass){
	RooAbsPdf* pdf = 0;
	for(int i=0; i<lines.size(); i++){
		if(c==lines[i][2] and p==lines[i][1]){
			RooWorkspace *w=0;
			w = (RooWorkspace*)GetRWSfromMap(MAP_RWSname_Pointer, lines[i][3], GetWordFromLine(lines[i][4],0,":").Data());
			if(w==0)
			w = (RooWorkspace*)GetTObject(lines[i][3], GetWordFromLine(lines[i][4], 0, ":").Data());
			AddRWSintoMap(lines[i][3], GetWordFromLine(lines[i][4],0,":").Data(), w, MAP_RWSname_Pointer);
			
			if(w==NULL) {
				cout<<" RooWorkspace is not found "<<endl; exit(1);
			}
			if(w->var("MH") && mass>0) w->var("MH")->setVal(mass);
			pdf= (RooAbsPdf*)w->pdf(GetWordFromLine(lines[i][4], 1 ,":")); //pdf->SetName(channelnames[c]+pdf->GetName());
			if(pdf==NULL) {
				w->Print();
				cout<<"pdf name = ["<<GetWordFromLine(lines[i][4], 1 ,":")<<"]"<<endl;
				cout<<" pdf is not found "<<endl; exit(1);
			}
		}
	}
	if(pdf==0) {
		for(int i=0; i<lines.size(); i++){
			for(int j=0; j<lines[i].size(); j++)cout<<lines[i][j]<<" ";
			cout<<""<<endl;
		}

		cout<<" Exist Workspaces : "<<endl;
		std::map<TString, RooWorkspace*>::iterator p1;
		for(p1=MAP_RWSname_Pointer.begin(); p1!=MAP_RWSname_Pointer.end(); ++p1){
			cout<<" in "<<p1->first<<endl;
		}
		cout<<"ERROR can't find TObject of process "<<p<<" in channel "<<c<<endl;
		exit(0);
	}
	return pdf;
}
RooAbsArg * GetExtraNorm(string c, string p, const vector< vector<string> > & lines, double mass){
	RooAbsArg * arg = 0;
	for(int i=0; i<lines.size(); i++){
		if(c==lines[i][2] and p==lines[i][1]){
			RooWorkspace *w=0;
			w = (RooWorkspace*)GetRWSfromMap(MAP_RWSname_Pointer, lines[i][3], GetWordFromLine(lines[i][4],0,":").Data());
			if(w==0)
			w = (RooWorkspace*)GetTObject(lines[i][3], GetWordFromLine(lines[i][4], 0, ":").Data());
			AddRWSintoMap(lines[i][3], GetWordFromLine(lines[i][4],0,":").Data(), w, MAP_RWSname_Pointer);
			if(w->var("MH") && mass>0) w->var("MH")->setVal(mass);
			arg= (RooAbsArg*)w->arg(GetWordFromLine(lines[i][4], 1 ,":")+"_norm"); //arg->SetName(channelnames[c]+arg->GetName());
		}
	}
	return arg;
}
RooAbsData* GetRooAbsData(string c, string p, const vector< vector<string> >& lines, double mass){
	RooAbsData* pdf = 0;
	for(int i=0; i<lines.size(); i++){
		if(c==lines[i][2] and p==lines[i][1]){
			RooWorkspace *w=0;
			w = (RooWorkspace*)GetRWSfromMap(MAP_RWSname_Pointer, lines[i][3], GetWordFromLine(lines[i][4],0,":").Data());
			if(w==0)
			w = (RooWorkspace*)GetTObject(lines[i][3], GetWordFromLine(lines[i][4], 0, ":").Data());
			AddRWSintoMap(lines[i][3], GetWordFromLine(lines[i][4],0,":").Data(), w, MAP_RWSname_Pointer);
			if(w->var("MH") && mass>0) w->var("MH")->setVal(mass);
			pdf= (RooAbsData*)w->data(GetWordFromLine(lines[i][4], 1 ,":")); //pdf->SetName(channelnames[c]+pdf->GetName());
		}
	}
	if(pdf==0) {
		cout<<"ERROR can't find TObject of process "<<p<<" in channel "<<c<<endl;
		exit(0);
	}
	return pdf;
}
void ReadLimitVsCLsFromFile(TGraphErrors*tge, TFile*f, int _debug) {
	tge->Set(0);	
	if (!f) {cout<<"Input file: "<<f->GetName()<<" error!  either not exist or empty,  exit"<<endl; exit(1);}
	TDirectory *toyDir = f->GetDirectory("");
	TIter next(toyDir->GetListOfKeys()); TKey *k;
	std::map<double, TTree*> gridCLsb; //r, <clsb, clserr>
	std::map<double, TTree*> gridCLb; //r, <clsb, clserr>
	std::map<double, double> gridQdata; //r, q_data
	while ((k = (TKey *) next()) != 0) {
		double rVal;
		TString name(k->GetName());
		if(name.BeginsWith("SAMPLING_SB_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_SB_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLsb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tsb->entries= "<<t->GetEntries()<< std::endl;
		}else if(name.BeginsWith("SAMPLING_B_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_B_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tb->entries= "<<t->GetEntries()<< std::endl;
		}else if(name.BeginsWith("DATA_R")){
			name.ReplaceAll("DATA_R","");
			TString tmp = name;
			rVal = tmp.Remove(tmp.Index("_"), tmp.Length()).Atof();
			name.Remove(0,name.Index("_Q")+2);
			gridQdata[rVal]=name.Atof();
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << tmp << " --> " << rVal << " Q_data ="<<name<<" --> "<<name.Atof()<< std::endl;
		}

	}

	int i = 0, n = gridCLsb.size();
	if(_debug)cout<<" grid size = "<<n<<endl;
	for (std::map<double, TTree *>::iterator itg = gridCLsb.begin(), edg = gridCLsb.end(); itg != edg; ++itg) {
		double cls, clserr;
		GetCLs(gridQdata[itg->first], itg->second, gridCLb[itg->first], cls, clserr);
		tge->SetPoint(     i, itg->first, cls   ); 
		tge->SetPointError(i, 0,          clserr);
		if(_debug>=10)cout<<" input grid:  r="<<itg->first<<" cls="<<cls<<"+/-"<<clserr<<endl;
		i++;
	}
	if(_debug)cout<<"tge->N = "<<tge->GetN()<<endl;
}
bool GetCLs(double qdata, TTree* tsb, TTree*tb,  double &cls, double &err, int _debug){
	vector<double> vclsb, vclb; vclsb.clear(); vclb.clear();
	GetM2lnQ(tsb, tb, vclsb, vclb);

	double clsb, clb, errsb, errb;
	GetPValue(vclsb, qdata, clsb, errsb);
	GetPValue(vclb, qdata, clb, errb);

	if(clb==0){if(_debug)	cout<<"CLsBase::CLs  Warning clb_b==0 !!!!"<<endl; err = 1; cls=1; return true;}
	err = sqrt( errb/clb*errb/clb + errsb/clsb*errsb/clsb) * clsb/clb;
	if(_debug>=10) cout<<"CLsBase::CLs  CLs=CLsb/CLb="<<clsb/clb<<"+/-"<<err<<endl;
	cls = clsb/clb;

	return true;
}
bool GetM2lnQ(TTree* tsb, TTree*tb, vector<double> &vclsb, vector<double>&vclb, int _debug){
	double clsb, clb;
	TBranch *brCLsb ;
	tsb->SetBranchAddress("brT", &clsb, &brCLsb);
	Long64_t nentries = tsb->GetEntries();
	for(int i=0; i<nentries; i++){
		tsb->GetEntry(i);
		vclsb.push_back(clsb);
	}
	//for(int i=0; i<vclsb.size(); i++) cout<<" "<<vclsb[i]<<" "<<endl;
	TBranch *brCLb ;
	tb->SetBranchAddress("brT", &clb, &brCLb);
	nentries = tb->GetEntries();
	for(int i=0; i<nentries; i++){
		tb->GetEntry(i);
		vclb.push_back(clb);
	}
	//for(int i=0; i<vclb.size(); i++) cout<<" "<<vclb[i]<<" "<<endl;
	return true;		
}
bool GetPValue(vector<double> vclsb, double qdata, double &ret, double &err, int _debug){
	double tmp = qdata;
	int _nexps = vclsb.size();
	for(int i=0; i<_nexps; i++){
		if(vclsb[i]>=tmp)
			ret ++ ;	
	}		
	ret/=_nexps;

	err= sqrt(ret*(1-ret)/_nexps);
	if(ret==0||ret==1) err= 1./_nexps;

	if(_debug>=10){
		cout<<"p ="<<ret<<" +/- "<<err<<" and total exps="<<_nexps<<endl;
		if(ret*_nexps <= 20) cout<<"p*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
	}
/*	if(ret == 0){
		if(_debug)	cout<<"p=0, it means number of pseudo experiments is not enough"<<endl;
		if(_debug)	cout<<"              Currently, we put p=1./"<<_nexps<<endl;
		ret = 1./(double)_nexps;
	}
*/
	return true;
}
void ReadM2lnQGridFromFile(TString filename, std::map<double, TTree*>&gridCLsb, std::map<double, TTree*>&gridCLb, int _debug) {
	if( gSystem->AccessPathName(filename)) {cout<<"file: ("<<filename<<") couldn't be found"<<endl; exit(0);};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
	if(f==NULL) f=new TFile(filename);
	if (!f) {cout<<"Input file: "<<f->GetName()<<" error!  either not exist or empty,  exit"<<endl; exit(1);}
	TDirectory *toyDir = f->GetDirectory("");
	TIter next(toyDir->GetListOfKeys()); TKey *k;
	while ((k = (TKey *) next()) != 0) {
		double rVal;
		TString name(k->GetName());
		if(name.BeginsWith("SAMPLING_SB_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_SB_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLsb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tsb->entries= "<<t->GetEntries()<< std::endl;
		}else if(name.BeginsWith("SAMPLING_B_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_B_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tb->entries= "<<t->GetEntries()<< std::endl;
		}
	}
	if(_debug) cout<<"M2lnQGrid size = "<<gridCLsb.size()<<endl;
	if(gridCLsb.size()!=gridCLb.size()) cout<<"Error: gridCLsb size="<<gridCLsb.size()<<" *** != *** "<<"gridCLb ="<<gridCLb.size()<<endl;
}

void ReadM2lnQGridFromFile(TString filename, std::map<double, double>&gridQdata, int _debug) {
	if( gSystem->AccessPathName(filename)) {cout<<"file: ("<<filename<<") couldn't be found"<<endl; exit(0);};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
	if(f==NULL) f=new TFile(filename);
	if (!f) {cout<<"Input file: "<<f->GetName()<<" error!  either not exist or empty,  exit"<<endl; exit(1);}
	TDirectory *toyDir = f->GetDirectory("");
	TIter next(toyDir->GetListOfKeys()); TKey *k;
	while ((k = (TKey *) next()) != 0) {
		double rVal;
		TString name(k->GetName());
		if(name.BeginsWith("DATA_R")){
			name.ReplaceAll("DATA_R","");
			TString tmp = name;
			rVal = tmp.Remove(tmp.Index("_"), tmp.Length()).Atof();
			name.Remove(0,name.Index("_Q")+2);
			gridQdata[rVal]=name.Atof();
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << tmp << " --> " << rVal << " Q_data ="<<name<<" --> "<<name.Atof()<< std::endl;
		}
	}
	if(_debug) cout<<"M2lnQGrid size = "<<gridQdata.size()<<endl;
}
vector<double> GetVectorFrom(TTree* tree, TString brName){
	vector<double> v;
	double c;
	TBranch *br;
	tree->SetBranchAddress(brName, &c, &br);
	Long64_t nentries = tree->GetEntries();
	for(int i=0; i<nentries; i++){
		tree->GetEntry(i);
		v.push_back(c);
	}
	return v;
}

RooWorkspace* GetRWSfromMap(map<TString,RooWorkspace*>m, string filename, string rwsname){
	RooWorkspace *w=0;
	std::map<TString, RooWorkspace*>::iterator p;
	TString s = "FILENAME_"; s+=filename; s+="_RWSNAME_"; s+=rwsname;
	//cout<<" To be found = "<<s<<endl;
	for(p=m.begin(); p!=m.end(); ++p){
		//cout<<" in "<<p->first<<endl;
		if(s==p->first) { 
			//cout<<" Exist RWS in memory, take it"<<endl;
			return p->second;	
		}
	}
			//cout<<" NOT Exist RWS in memory"<<endl;
	return w;
}
void AddRWSintoMap(string filename, string rwsname, RooWorkspace* w, map<TString,RooWorkspace*>& m ){
	TString s = "FILENAME_"; s+=filename; s+="_RWSNAME_"; s+=rwsname;
	m[s]=w;
}
