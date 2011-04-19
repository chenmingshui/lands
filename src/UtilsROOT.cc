#include "UtilsROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "RooDataSet.h"
TH1F* GetHisto(string filename, string histoname){

	//cout<<filename<<", "<<histoname<<endl;
	// FIXME need to check if filename is exist, and histoname is exist 
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; exit(0);};
	TFile *f =new TFile(filename.c_str());
	TH1F *h = (TH1F*)f->Get(histoname.c_str());
	if(!h) {cout<<"hist ["<<histoname<<"] in file ["<<filename<<"] couldn't be found"<<endl; exit(0);};
	return h;
}
TObject* GetTObject(string filename, string objname){

	// FIXME need to check if filename is exist, and histoname is exist 
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; exit(0);};
	TFile *f =new TFile(filename.c_str());
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
void FillTree(TString sfile, vector<double> array){
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
int GetTotProc(int c, int p, vector< vector<string> >vv_procnames){
	int stop = 0;
	for(int i=0; i<c; i++){
		stop += vv_procnames[i].size();
	}
	return stop+p;
}
string GetUncertainy(string c, vector<string> channelnames, string p, vector< vector<string> >vv_procnames, vector<string>ss1){
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
string GetUncertainy(int c, int p, vector< vector<string> >vv_procnames, vector<string>ss1){
	int stop = 0;
	for(int i=0; i<c; i++){
		stop += vv_procnames[i].size();
	}
	//cout<<stop+p+3<<endl;
	if(ss1[1]=="gmN" or ss1[1]=="gmA") return ss1[stop+p+3];
	else return ss1[stop+p+2];
}
bool CheckIfDoingShapeAnalysis(CountingModel* cms, TString ifileContentStripped, int debug){
	//int debug = 0;
	vector<TString> lines;
	lines = SplitIntoLines(ifileContentStripped, debug);

	// check if there is key word "shape", if yes, we need expand the
	// whole file to include all bins of shapes
	bool hasShape = false;
	bool hasParametricShape = false;
	vector<TString> shapeinfo; 
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("shapes")){
			hasShape = true;
			shapeinfo.push_back(lines[l]);
			if(GetWordFromLine(lines[l], 4).Contains(":")) hasParametricShape =true;
		}
	}
	if(hasShape==false) ConfigureModel(cms, ifileContentStripped, debug);
	if(hasShape){
		// get number of channels
		int nchannel = -1;
		for(int l=0; l<lines.size(); l++){
			if(lines[l].BeginsWith("imax")){
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
			if(lines[l].BeginsWith("Observation") or lines[l].BeginsWith("observation")){
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
			if(lines[l].BeginsWith("process")){
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
			if(lines[l].BeginsWith("rate")){
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
			if(lines[l].BeginsWith("kmax")){
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
		for(int j=lines.size()-1; j>=0; j--){
			if(lines[j].BeginsWith("---") )	break;
			if(lines[j].BeginsWith("rate"))	break;
			if(lines[j].BeginsWith("bin") )	break;
			if(lines[j].BeginsWith("process"))break;
			if(lines[j].BeginsWith("imax"))break;
			if(lines[j].BeginsWith("jmax"))break;
			if(lines[j].BeginsWith("kmax"))break;
			if(lines[j].BeginsWith("shapes"))continue;
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[j].Data(), " ");
			if(ss.size()<ntotprocesses+2 and GetWordFromLine(lines[j], 1)!="param" ){
				cout<<"uncertainty configuration is not correct"<<endl; 
				cout<<lines[j]<<endl;
				exit(0);
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
		for(int j=0; j<shapeinfo.size(); j++){
			if(shapeinfo[j].BeginsWith("shapes")){
				vector<string>	ss; 
				ss.clear();
				StringSplit(ss, shapeinfo[j].Data(), " ");
				if(ss.size()>=5){
					if(ss[1]=="*" || ss[2]=="*") xx_needtointerpret.push_back(ss);
					else {
						shape.push_back(ss);
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
							shape.push_back(newline);
						}
					}
				}
			}
		}
		for(int i=0; i<xx_needtointerpret.size(); i++){
			if(xx_needtointerpret[i][INDEXofChannel]=="*" or xx_needtointerpret[i][INDEXofProcess]=="*"){
				for(int c=0; c<nchannel; c++){
					if(xx_needtointerpret[i][INDEXofChannel]=="*" or xx_needtointerpret[i][INDEXofChannel]==channelnames[c]){
						for(int p=0; p<vv_procnames[c].size(); p++){
							if(xx_needtointerpret[i][INDEXofProcess]=="*" or xx_needtointerpret[i][INDEXofProcess]==vv_procnames[c][p]){
								vector<string> newline = xx_needtointerpret[i];
								newline[INDEXofChannel] = channelnames[c]; newline[INDEXofProcess]=vv_procnames[c][p];
								newline[4] = TString(newline[4]).ReplaceAll("$CHANNEL", channelnames[c]);
								newline[4] = TString(newline[4]).ReplaceAll("$PROCESS", vv_procnames[c][p]);
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
		}
		for(int i=0; i<shapeuncertainties.size(); i++){// for histogram morphing ....    keyword:  shapeN, shapeL, shape, shapeQ ,  alternative two sets of histogram for the variations w.r.t a particular source of systematics
			for(int j=i+1; j<shapeuncertainties.size(); j++){
				if(shapeuncertainties[j][1]==shapeuncertainties[i][1] and shapeuncertainties[j][2]==shapeuncertainties[i][2]
						and shapeuncertainties[j][5] == shapeuncertainties[i][5]
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


		int newchannels = nchannel;
		for(int c = 0; c<channelnames.size(); c++){
			for(int j=0; j<shape.size(); j++){
				TString s = shape[j][INDEXofChannel];
				if(s!=channelnames[c])continue;
				if(TString(shape[j][4]).Contains(":")) continue; // parametric channel, will process later
				TH1F *h = (TH1F*)GetHisto(shape[j][3],shape[j][4]);
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
			if(lines[k].BeginsWith("Observation") or lines[k].BeginsWith("observation")){
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
			if(ss1[1]=="shapeN") ss1[1]="lnN";
			TString s = ""; s+=ss1[0].c_str(); s+=" "; s+=ss1[1].c_str(); s+=" "; // weired, it doesn't give good error message when using "s+=ss1[1]"
			if(ss1[1]=="gmA" or ss1[1]=="gmN"){
				s+=ss1[2]; s+=" ";
				//cout<<s<<endl;
			};
			//if(uncerlines[u][1]=="gmN") cout<<s<<endl;
			if(uncertypes[u]=="param") {
				for(int ii = 2; ii<ss1.size(); ii++){
					s+=ss1[ii]; s+=" ";
				}
			}
			vs_unc.push_back(s);
		}

		int n = 0; // for index of new whole channels
		for(int c=0; c<nchannel; c++){//old channel number will be replaced 
			bool isShapeChannel = false;
			for(int p=0; p<shape.size(); p++){
				if(shape[p][INDEXofChannel]!=channelnames[c] || shape[p][INDEXofProcess]!="data_obs") continue; // find the signal one
				if(TString(shape[p][4]).Contains(":")) continue; // parametric channel, will process later

				isShapeChannel=true;

				if(debug)cout<<"channel="<<shape[p][INDEXofChannel]<<" data_obs file="<<shape[p][3]<<" hist="<<shape[p][4]<<endl;

				int n_proc = nprocesses[c];
				TH1F* h = (TH1F*)GetHisto(shape[p][3], shape[p][4]);
				TH1F* hn[n_proc];
				TH1F *hnorm[n_proc];
				TH1F *hunc_up[n_proc][nsyssources];
				TH1F *hunc_dn[n_proc][nsyssources];
				TH1F *hunc_up_norm[n_proc][nsyssources];
				TH1F *hunc_dn_norm[n_proc][nsyssources];
				for(int t=0; t<n_proc; t++){
					for(int q=0; q<shape.size(); q++){
						if(shape[q][INDEXofChannel]!=channelnames[c] || shape[q][INDEXofProcess]!=vv_procnames[c][t]) continue; 
						hn[t] = ((TH1F*)GetHisto(shape[q][3], shape[q][4]));
						if(hn[t]->Integral()== 0) { 
							cout<<"ERROR: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
								<<"] is shape, but the norminal histogram->Integral = 0"<<endl;
							exit(0);
						}
						hnorm[t]=(TH1F*)hn[t]->Clone("del_clone"+t);
						hnorm[t]->Scale(1./hn[t]->Integral());
					}

					// for channel c, process t,  looking for uncer which is shaping 
					for(int u=0; u<nsyssources; u++){
						TString type = uncerlines[u][1];
						if(type.BeginsWith("shape")){ 
							TString unc = GetUncertainy(c, t, vv_procnames, uncerlines[u]);
							if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
								bool filled = false;
								for(int i=0; i<shapeuncertainties.size(); i++){ // look for the histograms from the list 
									if(shapeuncertainties[i][INDEXofChannel]==channelnames[c] 
											and shapeuncertainties[i][INDEXofProcess]==vv_procnames[c][t]
											and shapeuncertainties[i][4]==uncerlines[u][0]){
										hunc_up[t][u] = (TH1F*)GetHisto(shapeuncertainties[i][3], shapeuncertainties[i][5]);
										if(hunc_up[t][u]->Integral()== 0) { 
											cout<<"ERROR: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
												<<"] is shape, but the up_shift histogram->Integral = 0"<<endl;
											exit(0);
										}
										hunc_dn[t][u] = (TH1F*)GetHisto(shapeuncertainties[i][3], shapeuncertainties[i][6]);
										if(hunc_dn[t][u]->Integral()== 0) { 
											cout<<"ERROR: channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]
												<<"] is shape, but the down_shift histogram->Integral = 0"<<endl;
											exit(0);
										}
										hunc_up_norm[t][u]=(TH1F*)hunc_up[t][u]->Clone("up_clone"+t+u);
										hunc_up_norm[t][u]->Scale(1./hunc_up[t][u]->Integral());
										hunc_dn_norm[t][u]=(TH1F*)hunc_dn[t][u]->Clone("dn_clone"+t+u);
										hunc_dn_norm[t][u]->Scale(1./hunc_dn[t][u]->Integral());
										filled = true;
									}
								}	
								if(!filled) {
									cout<<"ERROR channel ["<<channelnames[c]<<"] process ["<<vv_procnames[c][t]<<"] sys ["<<uncerlines[u][0]
										<<"] is shape uncertainty, we need corresponding histogram "<<endl;
									exit(0);
								}
							}
						}
					}
				}



				for(int r=1; r<=h->GetNbinsX(); r++){
					n++; 
					int proc = 0;
					for(int t=0; t<n_proc; t++){
						s_bin+=channelnames[c]; s_bin+="_"; s_bin+=r; s_bin+=" ";
						s_rate+=((hn[proc]))->GetBinContent(r);s_rate+=" ";
						// fill  s_process  //FIXME
						s_process += (proc-nsigproc[c]+1) ; s_process += " ";
						s_process_name += vv_procnames[c][proc] ; s_process_name += " ";

						if(debug>=10)cout<<"channel "<<c <<", histo bin"<<r<<" oldbin "<<t<<endl;

						for(int u=0; u<nsyssources; u++){
							//for(int k=0; k<uncerlines.size(); k++)
							if(uncertypes[u]=="param") continue;
							if(1){
								if(debug>=100) cout<<"debug "<<uncertypes[u]<<endl;
								if(uncertypes[u]=="shape" or uncertypes[u]=="shapeL" or uncertypes[u]=="shapeQ" or uncertypes[u]=="shapeN"){
									TString unc = GetUncertainy(c, t, vv_procnames, uncerlines[u]);
									if(unc.IsFloat() && unc.Atof()>0){ // number should be > 0
										double down = hunc_dn_norm[t][u]->GetBinContent(r);
										if(debug) cout<<"down "<<down<<endl;
										double up = hunc_up_norm[t][u]->GetBinContent(r);
										double norminal = hnorm[t]->GetBinContent(r);
										double norm_down = hunc_dn[t][u]->Integral();
										double norm_up = hunc_up[t][u]->Integral();
										double norm_norminal = hn[t]->Integral();
										TString stmp;
										if(uncertypes[u]=="shapeN"){
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
									}else{
										vs_unc[u]+=" - ";
									}
								}else{
									string tmp =GetUncertainy(c, t, vv_procnames, uncerlines[u]); 
									if(uncertypes[u]=="gmA" or uncertypes[u]=="gmN"){
										cout<<"gamma "<<tmp<<endl;
										if(TString(tmp).IsFloat()){
											float tmpd = TString(tmp).Atof()*hnorm[t]->GetBinContent(r);
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
				for(int t=0; t<n_proc; t++){
					for(int u=0; u<nsyssources; u++){
						/*
						   if(hunc_up[t][u]) delete hunc_up[t][u];
						   if(hunc_up_norm[t][u]) delete hunc_up_norm[t][u];
						   if(hunc_dn[t][u]) delete hunc_dn[t][u];
						   if(hunc_dn_norm[t][u]) delete hunc_dn_norm[t][u];
						   */
					}
				}

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
						if(uncertypes[u]=="param") continue;
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
		//exit(1);
		if(hasParametricShape) ConfigureShapeModel(cms, cardExpanded, parametricShapeLines,  debug);
		else ConfigureModel(cms, cardExpanded, debug);
	}// line begin with "shape"
}

bool ConfigureModel(CountingModel *cms, TString ifileContentStripped, int debug){
	// channel index starts from 1
	// systematics source index also starts from 1

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
		if(lines[l].BeginsWith("imax")){
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
		if(lines[l].BeginsWith("Observation") or lines[l].BeginsWith("observation")){
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
		if(lines[l].BeginsWith("process")){
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
		if(lines[l].BeginsWith("rate")){
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
		if(lines[l].BeginsWith("kmax")){
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
		cms->AddChannel(channelnames[c], sigbkgs, nsigproc[c]);
		cms->SetProcessNames(c, tmpprocn);
		cms->AddObservedData(c, observeddata[c]);
	}	

	//cms->Print();

	if(debug) cout<<"filled yields"<<endl;
	// fill the model with systemics

	// sections are separated by "---------"  in the data card
	// last section are all uncertainties,  it's analyzers' responsibility to make it right  
	int nsyssources = 0;
	for(int j=lines.size()-1; j>=0; j--){
		if(lines[j].BeginsWith("---") )	break;
		if(lines[j].BeginsWith("rate"))	break;
		if(lines[j].BeginsWith("bin") )	break;
		if(lines[j].BeginsWith("process"))break;
		if(lines[j].BeginsWith("imax"))break;
		if(lines[j].BeginsWith("jmax"))break;
		if(lines[j].BeginsWith("kmax"))break;
		nsyssources++;
	}
	if(kmax>=0 && kmax!=nsyssources) {cout<<"kmax !=  number of independant uncertainties"<<endl; exit(0);}

	for(int s=0; s<nsyssources; s++){
		TString tmps = TString::Format("%6d ",s+1);
		vector<string>	ss; 

		int hasUncSourceWithIndex_s = -1;
		/*
		   for(int j=0; j<lines.size(); j++){
		   ss.clear();
		   StringSplit(ss, lines[j].Data(), " ");
		   if (TString(ss[0]).Atoi()==(s+1)) {
		   if(hasUncSourceWithIndex_s >= 0) cout<<"Warning:  There are more than one uncertainty source start with index = "<<s+1<<", will use last one"<<endl; 
		   hasUncSourceWithIndex_s = j;
		   }
		   }
		   if(hasUncSourceWithIndex_s<0) {
		   cout<<"Error: no uncertainty start with index = "<<s+1<<endl;
		   cout<<"   Prob. reasons: check kmax or index of systematics sources"<<endl;
		   exit(0); 
		   }
		   */
		hasUncSourceWithIndex_s = lines.size() - nsyssources + s;
		ss.clear();
		StringSplit(ss, lines[hasUncSourceWithIndex_s].Data(), " ");

		/*
		   if(TString(ss[0]).IsDigit() == false) {
		   cout<<"   You are trying to assign the following line which not starting from digits to "<< s+1 <<"th systematics sources"<<endl;
		   cout<<"\""<<lines[hasUncSourceWithIndex_s]<<"\""<<endl;
		   cout<<"   Probable reasons: check kmax or index of systematics sources"<<endl;
		   exit(0);
		   }
		   */
		// FIXME need to check: don't have two lines with same uncertainty name...

		//int indexcorrelation = (TString(ss[0])).Atoi();
		string indexcorrelation = ss[0];

		//if(indexcorrelation != (nsyssources-s) ) {
		//	cout<<"Warning: you are reading in "<<nsyssources-s<<"th nuisanse parameter, but this line start with "<<indexcorrelation<<endl;
		//}

		int pdf = 0; 
		// see CountingModel.h
		if(ss[1]=="lnN") pdf=typeLogNormal;   // typeLogNormal
		else if(ss[1]=="trG") pdf=typeTruncatedGaussian; // typeTruncatedGaussian
		else if(ss[1]=="gmA" or ss[1]=="gmN" or ss[1]=="gmM") pdf=typeGamma; // typeControlSampleInferredLogNormal;  gmA was chosen randomly, while in gmN, N stands for yield in control sample;  gmM stands for Multiplicative gamma distribution, it implies using a Gamma distribution not for a yield but for a multiplicative correction
		else if(ss[1]=="shapeN") pdf=typeLogNormal;
		else if(ss[1]=="shape" or ss[1]=="shapeQ") pdf=typeShapeGaussianQuadraticMorph;
		else if(ss[1]=="shapeL" ) pdf=typeShapeGaussianLinearMorph;
		else pdf =  (TString(ss[1])).Atoi();

		tmps+= TString::Format("%3s ", ss[1].c_str());
		if(pdf==typeGamma && ss[1]!="gmM") tmps+= TString::Format("%8s ", ss[2].c_str());
		else		   tmps+= "         ";

		bool filledThisSource = false;
		for(int p=0; p<ntotprocesses; p++){
			double err, errup, rho;  // to allow asymetric uncertainties
			double shape[8];
			if(pdf==typeLogNormal){
				if(ss[p+2]=="-") {
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/")){
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, ss[p+2], "/");
					if((TString(asymetricerrors[0])).Atof()<=0) { cout<<"ERROR:  Kappa can't be <=0 "<<endl; exit(0); };
					err= 1./(TString(asymetricerrors[0])).Atof()-1.0; // downside 
					errup= (TString(asymetricerrors[1])).Atof()-1.0;  // upside
				}else {
					err= (TString(ss[p+2])).Atof()-1.0;
					errup = err;
				}
				tmps+= TString::Format("%7.2f ",err+1.0);
				if(err < -1 or errup <-1) {
					cout<<"ERROR: Kappa in lognormal  can't be negative, please check your input card at "<<s+1<<"th source, "<<p+2<<"th entry"<<endl;
					exit(0);
				}
				if(err == 0. && errup==0.) continue;
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
				if(TString(ss[p+2]).Contains("/")){
					vector<string> params; params.clear();
					StringSplit(params, ss[p+2], "/");
					if(params.size()==7) {
						for(int i=0; i<7; i++) shape[i]=TString(params[i]).Atof();
						shape[7]=1.;
					}else if(params.size()==2){
						for(int i=0; i<6; i++) {
							if(i<4)shape[i]= 0 ;
							else shape[i]=TString(params[i-4]).Atof();
						}
						shape[6]=0;
						shape[7]=0;
					}else { cout<<"ERROR: typeShape input format not correct "<<endl; exit(0);}
				}else{
					for(int i=0; i<4; i++) {
						shape[i]= 0 ;
					}
					shape[4]=TString(ss[p+2]).Atof();
					shape[5]=shape[4];
					shape[6]=0;
					shape[7]=0;
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
			//cout<<"delete me: c="<<binnumber[p]-1<<" s="<< subprocess[p]<<endl;
			if(pdf==typeLogNormal||pdf==typeTruncatedGaussian)cms->AddUncertainty(binnumber[p]-1, subprocess[p], err, errup, pdf, indexcorrelation );
			if(pdf==typeGamma){
				if(ss[1]=="gmA" or ss[1]=="gmN"){
					double N = (TString(ss[2]).Atof()); // your input should be Most Probable Value,  while mean value is N+1
					if(N<0)  {
						cout<<"Yield in control Region in gamma can't be negative, please check your input card at "<<s+1<<"th source, "<<2<<"th entry"<<endl;
						exit(0);
					}
					cms->AddUncertainty(binnumber[p]-1, subprocess[p], rho, 0, N, pdf, indexcorrelation );
					//FIXME  here need to check:  rho*N == cms->GetExpectedNumber(binnumber[p]-1, subprocess[p]);
				}
				else if(ss[1]=="gmM"){
					double N = 1./err/err-1; // we use convention: mean of events is 1./err/err,  so we do "-1" here, and then add back "1"  in src/CountingModel.cc
					//double N = (int) (1./err/err);
					cms->AddUncertainty(binnumber[p]-1, subprocess[p], -1, 0, N, pdf, indexcorrelation ); // use "rho = -1" here to imply that this gamma term is not for control sample inferred uncertainties, but multiplicative gamma function ....
					//FIXME  here need to check:  all uncertainties with same indexcorrelation are the same,   can't be different currently ... 
					//  e.g.   check   N == cms->Get_v_Gamma(indexcorrelation); // can't do here before ConfigUncertaintyPdfs()
					//we might allow them different and do rescaling 
				}
			}
			if(pdf==typeShapeGaussianLinearMorph or pdf==typeShapeGaussianQuadraticMorph){
				cms->AddUncertainty(binnumber[p]-1, subprocess[p], 8, shape, pdf, indexcorrelation );
			}

			// because when err < 0, AddUncertainty do nothing,  but filledThisSource has been changed to be true
			filledThisSource = true;
		}
		if(!filledThisSource) {
			cout<<"WARNING: The "<< s+1 <<"th source of uncertainties are all 0. "<<endl;
			//FIXME temporarily solution:  when reading a source with all error = 0,  then assign it to be logNormal, error =0,  in UtilsROOT.cc 
			// to incorporate the MinuitFit.  This ad-hoc fix will slow down the toy generation because there are lots of non-use random numebrs and operations ...
			// cms->AddUncertainty(0, 0, 0, 1, indexcorrelation );
			//cms->AddUncertainty(0, 0, 0, 1, indexcorrelation ); //FIXME  no need now, because we use name of uncertainties, 
			//exit(0);
		}
		duplicatingLines.push_back(tmps);
	}

	if(debug) cout<<"filled systematics"<<endl;

	//cms->Print(100);
	if(debug) {
		cout<<"start duplicating this card:"<<endl<<endl;
		for(int i=0; i<duplicatingLines.size(); i++) cout<<duplicatingLines[i]<<endl;
		cout<<endl<<endl<<"end duplicating this card:"<<endl<<endl;
	}

	return true;
}

bool ConfigureModel(CountingModel *cms, const char* fileName, int debug){
	TString s = ReadFile(fileName);
	cms->SetModelName(fileName);
	return CheckIfDoingShapeAnalysis(cms, s, debug);
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
bool ConfigureShapeModel(CountingModel *cms, TString ifileContentStripped, vector< vector<string> > parametricShapeLines, int debug){
	// channel index starts from 1
	// systematics source index also starts from 1
	if(debug)for(int i=0; i<parametricShapeLines.size(); i++){
		for(int j=0; j<parametricShapeLines[i].size(); j++){
			cout<<parametricShapeLines[i][j]<<" ";
		}
		cout<<endl;
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
		if(lines[l].BeginsWith("imax")){
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
		if(lines[l].BeginsWith("Observation") or lines[l].BeginsWith("observation")){
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
		if(lines[l].BeginsWith("process")){
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
		if(lines[l].BeginsWith("rate")){
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
		if(lines[l].BeginsWith("kmax")){
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
		for(int i=0; i<parametricShapeLines.size(); i++){
			if(parametricShapeLines[i][1]==channelnames[c]);
			isParametricChannel = true;
		}
		if(!isParametricChannel){
			cms->AddChannel(channelnames[c], sigbkgs, nsigproc[c]);
			cms->SetProcessNames(c, tmpprocn);
			cms->AddObservedData(c, observeddata[c]);
		}else{
			vector<RooAbsPdf*> vspdf, vbpdf; vspdf.clear(); vbpdf.clear();
			for(int i=0; i<sigbkgs.size(); i++){
				RooAbsPdf * pdf = (RooAbsPdf*)GetPdf(channelnames[c], tmpprocn[i], parametricShapeLines);
				pdf->SetName(TString(channelnames[c])+pdf->GetName());
				if(i<nsigproc[c])vspdf.push_back(pdf);
				else vbpdf.push_back(pdf);
			}
			vector<double> vsnorm, vbnorm; vsnorm.clear(); vbnorm.clear();
			for(int i=0; i<sigbkgs.size(); i++){
				if(i<nsigproc[c]) vsnorm.push_back(sigbkgs[i]);
				else vbnorm.push_back(sigbkgs[i]);
			}
			RooDataSet *data = (RooDataSet*)GetRooDataSet(channelnames[c], "data_obs", parametricShapeLines);
			RooRealVar* x = dynamic_cast<RooRealVar*> (data->get()->first());
			x->SetName(TString(channelnames[c])+x->GetName());
			RooWorkspace *w;
			cms->AddChannel(channelnames[c], x, vspdf, vsnorm, vbpdf, vbnorm, w);
			cms->AddObservedDataSet(channelnames[c], data);
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
		if(lines[j].BeginsWith("rate"))	break;
		if(lines[j].BeginsWith("bin") )	break;
		if(lines[j].BeginsWith("process"))break;
		if(lines[j].BeginsWith("imax"))break;
		if(lines[j].BeginsWith("jmax"))break;
		if(lines[j].BeginsWith("kmax"))break;
		nsyssources++;
	}
	if(kmax>=0 && kmax!=nsyssources) {cout<<"kmax !=  number of independant uncertainties"<<endl; exit(0);}

	vector< vector<string> > vshape_params_unclines; vshape_params_unclines.clear();
	cout<<" total "<<nsyssources<<" systs"<<endl;
	for(int s=0; s<nsyssources; s++){
		vector<string>	ss; 

		int hasUncSourceWithIndex_s = -1;
		hasUncSourceWithIndex_s = lines.size() - nsyssources + s;
		ss.clear();
		StringSplit(ss, lines[hasUncSourceWithIndex_s].Data(), " ");
		if(debug)cout<<" @@@@@@@ systemticline: "<<lines[hasUncSourceWithIndex_s]<<endl;

		TString tmps; tmps.Form("%s ", ss[0].c_str());

		string indexcorrelation = ss[0];


		int pdf = 0; 
		// see CountingModel.h
		if(ss[1]=="lnN") pdf=typeLogNormal;   // typeLogNormal
		else if(ss[1]=="trG") pdf=typeTruncatedGaussian; // typeTruncatedGaussian
		else if(ss[1]=="gmA" or ss[1]=="gmN" or ss[1]=="gmM") pdf=typeGamma; // typeControlSampleInferredLogNormal;  gmA was chosen randomly, while in gmN, N stands for yield in control sample;  gmM stands for Multiplicative gamma distribution, it implies using a Gamma distribution not for a yield but for a multiplicative correction
		else if(ss[1]=="shapeN") pdf=typeLogNormal;
		else if(ss[1]=="shape" or ss[1]=="shapeQ") pdf=typeShapeGaussianQuadraticMorph;
		else if(ss[1]=="shapeL" ) pdf=typeShapeGaussianLinearMorph;
		else if(ss[1]=="param" ) {
			pdf=typeBifurcatedGaussian;
			vshape_params_unclines.push_back(ss);
		}
		else pdf =  (TString(ss[1])).Atoi();

		tmps+= TString::Format("%8s ", ss[1].c_str());
		if(pdf==typeGamma && ss[1]!="gmM") tmps+= TString::Format("%8s ", ss[2].c_str());
		else		   tmps+= "         ";

		bool filledThisSource = false;
		for(int p=0; p<ntotprocesses; p++){
			double err, errup, rho;  // to allow asymetric uncertainties
			double shape[8];
			if(pdf==typeLogNormal){
				if(ss[p+2]=="-") {
					tmps+= "    -   ";
					continue;
				}
				if(TString(ss[p+2]).Contains("/")){
					vector<string> asymetricerrors; asymetricerrors.clear();
					StringSplit(asymetricerrors, ss[p+2], "/");
					if((TString(asymetricerrors[0])).Atof()<=0) { cout<<"ERROR:  Kappa can't be <=0 "<<endl; exit(0); };
					err= 1./(TString(asymetricerrors[0])).Atof()-1.0; // downside 
					errup= (TString(asymetricerrors[1])).Atof()-1.0;  // upside
				}else {
					err= (TString(ss[p+2])).Atof()-1.0;
					errup = err;
				}
				tmps+= TString::Format("%7.2f ",err+1.0);
				if(err < -1 or errup <-1) {
					cout<<"ERROR: Kappa in lognormal  can't be negative, please check your input card at "<<s+1<<"th source, "<<p+2<<"th entry"<<endl;
					exit(0);
				}
				if(err == 0. && errup==0.) continue;
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
				if(TString(ss[p+2]).Contains("/")){
					vector<string> params; params.clear();
					StringSplit(params, ss[p+2], "/");
					if(params.size()==7) {
						for(int i=0; i<7; i++) shape[i]=TString(params[i]).Atof();
						shape[7]=1.;
					}else if(params.size()==2){
						for(int i=0; i<6; i++) {
							if(i<4)shape[i]= 0 ;
							else shape[i]=TString(params[i-4]).Atof();
						}
						shape[6]=0;
						shape[7]=0;
					}else { cout<<"ERROR: typeShape input format not correct "<<endl; exit(0);}
				}else{
					for(int i=0; i<4; i++) {
						shape[i]= 0 ;
					}
					shape[4]=TString(ss[p+2]).Atof();
					shape[5]=shape[4];
					shape[6]=0;
					shape[7]=0;
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
			//cout<<"delete me: c="<<binnumber[p]-1<<" s="<< subprocess[p]<<endl;
			if(pdf==typeLogNormal||pdf==typeTruncatedGaussian){
				cms->AddUncertaintyOnShapeNorm(binnumber[p]-1, subprocess[p], err, errup, pdf, indexcorrelation );
			}
			if(pdf==typeGamma){
				if(ss[1]=="gmA" or ss[1]=="gmN"){
					double N = (TString(ss[2]).Atof()); // your input should be Most Probable Value,  while mean value is N+1
					if(N<0)  {
						cout<<"Yield in control Region in gamma can't be negative, please check your input card at "<<s+1<<"th source, "<<2<<"th entry"<<endl;
						exit(0);
					}
					cms->AddUncertainty(binnumber[p]-1, subprocess[p], rho, 0, N, pdf, indexcorrelation );
					//FIXME  here need to check:  rho*N == cms->GetExpectedNumber(binnumber[p]-1, subprocess[p]);
				}
				else if(ss[1]=="gmM"){
					double N = 1./err/err-1; // we use convention: mean of events is 1./err/err,  so we do "-1" here, and then add back "1"  in src/CountingModel.cc
					//double N = (int) (1./err/err);
					cms->AddUncertainty(binnumber[p]-1, subprocess[p], -1, 0, N, pdf, indexcorrelation ); // use "rho = -1" here to imply that this gamma term is not for control sample inferred uncertainties, but multiplicative gamma function ....
					//FIXME  here need to check:  all uncertainties with same indexcorrelation are the same,   can't be different currently ... 
					//  e.g.   check   N == cms->Get_v_Gamma(indexcorrelation); // can't do here before ConfigUncertaintyPdfs()
					//we might allow them different and do rescaling 
				}
			}
			if(pdf==typeShapeGaussianLinearMorph or pdf==typeShapeGaussianQuadraticMorph){
				cms->AddUncertainty(binnumber[p]-1, subprocess[p], 8, shape, pdf, indexcorrelation );
			}

			// because when err < 0, AddUncertainty do nothing,  but filledThisSource has been changed to be true
			filledThisSource = true;
		}
		if(!filledThisSource) {
			cout<<"WARNING: The "<< s+1 <<"th source of uncertainties are all 0. "<<endl;
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
			cms->AddUncertaintyOnShapeParam(ss[0], mean, sigmaL, sigmaR, rangeMin, rangeMax );
		}

		duplicatingLines.push_back(tmps);
	}

	if(debug) cout<<"filled systematics"<<endl;

	//cms->Print(100);
	if(debug) {
		cout<<"start duplicating this card:"<<endl<<endl;
		for(int i=0; i<duplicatingLines.size(); i++) cout<<duplicatingLines[i]<<endl;
		cout<<endl<<endl<<"end duplicating this card:"<<endl<<endl;
	}

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
RooAbsPdf* GetPdf(string c, string p, vector< vector<string> > lines){
	RooAbsPdf* pdf = 0;
	for(int i=0; i<lines.size(); i++){
		if(c==lines[i][2] and p==lines[i][1]){
			RooWorkspace *w;
			w = (RooWorkspace*)GetTObject(lines[i][3], GetWordFromLine(lines[i][4], 0, ":").Data());
			pdf= (RooAbsPdf*)w->pdf(GetWordFromLine(lines[i][4], 1 ,":")); //pdf->SetName(channelnames[c]+pdf->GetName());
		}
	}
	if(pdf==0) {
		cout<<"ERROR can't find TObject of process "<<p<<" in channel "<<c<<endl;
		exit(0);
	}
	return pdf;
}
RooDataSet* GetRooDataSet(string c, string p, vector< vector<string> > lines){
	RooDataSet* pdf = 0;
	for(int i=0; i<lines.size(); i++){
		if(c==lines[i][2] and p==lines[i][1]){
			RooWorkspace *w;
			w = (RooWorkspace*)GetTObject(lines[i][3], GetWordFromLine(lines[i][4], 0, ":").Data());
			pdf= (RooDataSet*)w->data(GetWordFromLine(lines[i][4], 1 ,":")); //pdf->SetName(channelnames[c]+pdf->GetName());
		}
	}
	if(pdf==0) {
		cout<<"ERROR can't find TObject of process "<<p<<" in channel "<<c<<endl;
		exit(0);
	}
	return pdf;
}
