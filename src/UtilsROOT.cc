#include "UtilsROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
TH1F* GetHisto(string filename, string histoname){

	//cout<<filename<<", "<<histoname<<endl;

	// FIXME need to check if filename is exist, and histoname is exist 
	TFile *f =new TFile(filename.c_str());
	TH1F *h = (TH1F*)f->Get(histoname.c_str());
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
bool CheckIfDoingShapeAnalysis(CountingModel* cms, TString ifileContentStripped){
	int debug = 0;
	vector<TString> lines;
	lines = SplitIntoLines(ifileContentStripped, debug);
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
			if(s.IsDigit()==false) {
				cout<<"need a number after imax"<<", currently it's not a number but \'"<<s<<"\'"<<endl;
				exit(0);
			}
			nchannel = s.Atoi();
		} 
	}	
	if(nchannel==0) {
		cout<<"number of channels = 0,  please check imax"<<endl;
		return false;
	}

	//  get number of independant systematics sources
	int nsyssources = 0;
	bool hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("kmax")){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			nsyssources = (TString(ss[1])).Atoi();
			hasFilled = true;
		}
	}
	if (hasFilled == false) cout<<"need \"kmax\" to set number of independant systematic sources "<<endl;
	if (nsyssources==0) cout<<"no systematics input, kmax=0"<<endl;
	if(debug) cout<<"number of independant systematics sources = "<<nsyssources<<endl;

	// check if there is key word "shape", if yes, we need expand the
	// whole file to include all bins of shapes
	bool hasShape = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("shape")){
			hasShape = true;
			if(debug) cout<<"*************shape analysis**************"<<endl;
			TString cardExpanded;
			vector<TString> newlines; // will modify old file line by line
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");

			// must eliminate the "shape" in this new card,
			// otherwise it will enter into unlimited loop
			vector<TString> shapeinfo = SplitIntoLines(ReadFile(ss[1].c_str())); // additional file for shape infos

			if(debug) cout<<"shapeinfo lines = "<<shapeinfo.size()<<endl;

			vector<int> s_bins, s_procs, u_bins, u_procs, u_indexs;
			vector<TString> s_files, s_histos, u_files, u_histos;

			vector< vector<string> > shape;
			vector< vector<string> > uncertainty;
			vector< vector<string> > observation;
			for(int j=0; j<shapeinfo.size(); j++){
				if(shapeinfo[j].BeginsWith("shape")){
					vector<string>	ss; 
					ss.clear();
					StringSplit(ss, shapeinfo[j].Data(), " ");
					shape.push_back(ss);
				}
				if(shapeinfo[j].BeginsWith("uncertainty")){
					vector<string>	ss; 
					ss.clear();
					StringSplit(ss, shapeinfo[j].Data(), " ");
					uncertainty.push_back(ss);
				}
				if(shapeinfo[j].BeginsWith("observation")){
					vector<string>	ss; 
					ss.clear();
					StringSplit(ss, shapeinfo[j].Data(), " ");
					observation.push_back(ss);
				}
			}
			int newchannels = nchannel;
			for(int j=0; j<shape.size(); j++){
				TString s = shape[j][2];
				if(s.Atoi()!=0)continue;
				TH1F *h = (TH1F*)GetHisto(shape[j][3],shape[j][4]);
				newchannels+=h->GetNbinsX();
				newchannels-=1;
				delete h;
			}

			if(debug) cout<<"new number of channels: "<<newchannels<<endl;

			if(1){
				TString s = "imax "; s+=newchannels;
				newlines.push_back(s);
			}
			vector<string> ss_old_bin, ss_old_rate; 
			ss_old_bin.clear(); ss_old_rate.clear();
			for(int k=0; k<lines.size(); k++){
				if(lines[k].BeginsWith("kmax")){
					newlines.push_back(lines[k]);
				}
				if(lines[k].BeginsWith("Observation")){
					vector<string> ss_old; ss_old.clear();
					StringSplit(ss_old, lines[k].Data(), " ");
					TString s = "Observation ";
					int n = 0; // for index of new whole channels
					for(int c=1; c<=nchannel; c++){//old channel number will be replaced 
						bool isShapeChannel = false;
						for(int p=0; p<shape.size(); p++){
							if(TString(shape[p][1]).Atoi()!=c || TString(shape[p][2]).Atoi()!=0) continue; // find the signal one
							isShapeChannel=true;
							for(int q=0; q<observation.size(); q++){
								if(TString(observation[q][1]).Atoi()!=c) continue; // find the signal one
								TH1F* h = (TH1F*)GetHisto(observation[q][2], observation[q][3]);
								for(int r=1; r<=h->GetNbinsX(); r++){
									n++; s+=h->GetBinContent(r); s+=" ";
								}
								delete h;
							}
						}
						if(isShapeChannel==false) {n++; s+=ss_old[c]; s+=" ";}
					}	
					newlines.push_back(s);
				}
				if(lines[k].BeginsWith("bin")){
					StringSplit(ss_old_bin, lines[k].Data(), " ");
				}
				if(lines[k].BeginsWith("rate")){
					StringSplit(ss_old_rate, lines[k].Data(), " ");
				}
			}

			if(debug) cout<<"refill shape with observation"<<endl;

			TString s_bin = "bin ";
			TString s_rate = "rate ";
			vector<TString> vs_unc; vs_unc.clear();
			for(int u=1; u<=nsyssources; u++) {
				for(int k=0; k<lines.size(); k++){
					vector<string>	ss1; 
					ss1.clear();
					StringSplit(ss1, lines[k].Data(), " ");
					if(TString(ss1[0]).Atoi()!=u) continue;
					TString s = ""; s+=u; s+=" "; s+=ss1[1].c_str(); s+=" "; // weired, it doesn't give good error message when using "s+=ss1[1]"
					vs_unc.push_back(s);
				}
			}

			int n = 0; // for index of new whole channels
			for(int c=1; c<=nchannel; c++){//old channel number will be replaced 
				bool isShapeChannel = false;
				for(int p=0; p<shape.size(); p++){
					if(TString(shape[p][1]).Atoi()!=c || TString(shape[p][2]).Atoi()!=0) continue; // find the signal one
					isShapeChannel=true;

					if(debug)cout<<"channel "<<c <<", shape p "<<p<<endl;

					int n_proc = 0;
					for(int t=0; t<ss_old_bin.size(); t++){
						if(TString(ss_old_bin[t]).Atoi()!=c) continue;
						n_proc++;
					}
					TH1F* h = (TH1F*)GetHisto(shape[p][3], shape[p][4]);
					TH1F* hn[n_proc];
					//	= new TH1F[n_proc];
					for(int t=0; t<n_proc; t++){
						for(int q=0; q<shape.size(); q++){
							if(TString(shape[q][1]).Atoi()!=c || TString(shape[q][2]).Atoi()!=t) continue; // find the signal one
							hn[t] = ((TH1F*)GetHisto(shape[q][3], shape[q][4]));
						}
					}
					for(int r=1; r<=h->GetNbinsX(); r++){
						n++; 
						int proc = 0;
						for(int t=0; t<ss_old_bin.size(); t++){
							if(TString(ss_old_bin[t]).Atoi()!=c) continue;
							s_bin+=n; s_bin+=" ";
							s_rate+=((hn[proc]))->GetBinContent(r);s_rate+=" ";

							if(debug)cout<<"channel "<<c <<", histo bin"<<r<<" oldbin "<<t<<endl;

							for(int u=1; u<=nsyssources; u++){
								bool shapeuncertainty = false;
								for(int a=0; a<uncertainty.size(); a++){
									if(TString(uncertainty[a][1]).Atoi()==c && TString(uncertainty[a][2]).Atoi()==proc && TString(uncertainty[a][3]).Atoi()==u){
										shapeuncertainty = true;
										TH1F *hu = (TH1F*) GetHisto(uncertainty[a][4], uncertainty[a][5]);
										vs_unc[u-1]+=hu->GetBinContent(r); vs_unc[u-1] += " ";
										delete hu;
									}
								}
								if(shapeuncertainty==false){
									for(int k=0; k<lines.size(); k++){
										vector<string>	ss1; 
										ss1.clear();
										StringSplit(ss1, lines[k].Data(), " ");
										if(TString(ss1[0]).Atoi()!=u) continue;

										int totproc = 0;
										int proc_in_c = 0;
										for(int b=1; b<ss_old_bin.size(); b++){
											totproc++;
											if(TString(ss_old_bin[b]).Atoi()==c){
												if(proc_in_c==proc) break;
												proc_in_c++;
											}
										}
										vs_unc[u-1]+=ss1[totproc+1]; vs_unc[u-1] += " ";
									}
								}
							}



							proc++;
						}
					}
					delete h;
					for(int t=0; t<n_proc; t++) delete hn[t];

				}// shape channels
				if(isShapeChannel==false) {
					n++;
					int totproc = 0;
					for(int t=0; t<ss_old_bin.size(); t++){
						totproc++;
						if(TString(ss_old_bin[t]).Atoi()!=c) continue;
						s_bin+=n; s_bin+=" ";
						s_rate+=ss_old_rate[totproc];s_rate+=" ";
					}
				}
			}
			newlines.push_back(s_bin);
			newlines.push_back(s_rate);

			// extend uncertainty lines
			for(int u=1; u<=nsyssources; u++){
				// line start with int number  "u"		
				newlines.push_back(vs_unc[u-1]);
			}

			for(int j=0; j<newlines.size(); j++){
				cardExpanded+=newlines[j]; 
				if(j<(newlines.size()-1)) cardExpanded+="\n";
			}
			if(debug) cout<<"***********print out of the new huge table****************"<<endl;
			cout<<cardExpanded<<endl;
			//exit(1);
			ConfigureModel(cms, cardExpanded);

		}// line begin with "shape"
	}
	if(hasShape==false) ConfigureModel(cms, ifileContentStripped);
}

bool ConfigureModel(CountingModel *cms, TString ifileContentStripped){
	// channel index starts from 1
	// systematics source index also starts from 1

	TString duplicatingCard;
	vector<TString> duplicatingLines;

	bool debug = true;
	// Now proceed with the parsing of the stripped file

	vector<TString> lines;
	lines = SplitIntoLines(ifileContentStripped, debug);
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
			if(s.IsDigit()==false) {
				cout<<"need a number after imax"<<", currently it's not a number but \'"<<s<<"\'"<<endl;
				exit(0);
			}
			nchannel = s.Atoi();

			duplicatingLines.push_back(lines[l]);
		} 
	}	
	if(nchannel==0) {
		cout<<"number of channels = 0,  please check imax"<<endl;
		return false;
	}


	// get numbers of processes in each channel
	int *nprocesses = new int[nchannel];
	for(int i=0; i<nchannel; i++) nprocesses[i]=0;
	bool hasLineStartWithBin = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("bin") && hasLineStartWithBin==true){
			cout<<"WARNING:  there are two lines beginning with \"bin\" in your card. We will take the first one"<<endl;
		}
		if(lines[l].BeginsWith("bin") && hasLineStartWithBin==false){
			TString tmps = TString::Format("%10s ","bin");
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			for(int i=1; i<ss.size(); i++){
				int bin = (TString(ss[i])).Atoi();
				if(bin<(nchannel+1) && bin>=1) nprocesses[bin-1]++;
				else {
					cout<<"WARNING: there is a bin number = "<<bin<<", out of [1,"<<nchannel<<"]"<<endl;
					cout<<"We take "<<nchannel<<" bins indicated by imax, and skip the rest bins "<<endl;
					break;
					//return false;
				}
				tmps+= TString::Format("%7d ",bin);
			}
			hasLineStartWithBin = true;
			duplicatingLines.push_back(tmps);
		} 
	}	
	if(!hasLineStartWithBin) {
		cout<<"Line beginning with \"bin\" is not found in your card"<<endl;
		exit(0);
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

	// get observed dataset
	//
	double *observeddata=new double[nchannel];
	for(int i=0; i<nchannel; i++) observeddata[i]=0;
	hasFilled = false;
	for(int l=0; l<lines.size(); l++){
		if(lines[l].BeginsWith("Observation")){
			TString tmps = TString::Format("%10s ","Observation");
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
			if(ss.size()-1<nchannel) {cout<<"channels of observed data < "<<nchannel<<endl; return false;} 
			for(int i=1; i<(nchannel+1); i++){
				double ev = (TString(ss[i])).Atof();
				observeddata[i-1] = ev;
				tmps+= TString::Format("%7.2f ",ev);
			}
			hasFilled =  true;
			duplicatingLines.push_back(tmps);
		}
	}
	if(hasFilled==false) {cout<<"need a line starting with \"Observation\" "<<endl; return false;}
	if(debug){
		cout<<"observed data: ";
		for(int i=0; i<nchannel; i++) cout<<observeddata[i]<<" ";
		cout<<endl;
	}


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
			duplicatingLines.push_back(lines[l]);
		}
	}
	if (hasFilled == false) cout<<"need \"kmax\" to set number of independant systematic sources "<<endl;
	if (nsyssources==0) cout<<"no systematics input, kmax=0"<<endl;

	if(debug) cout<<"number of independant systematics sources = "<<nsyssources<<endl;


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
	// fill the model with systemics
	for(int s=0; s<nsyssources; s++){
		TString tmps = TString::Format("%6d ",s+1);
		vector<string>	ss; 

		int hasUncSourceWithIndex_s = -1;
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
		ss.clear();
		StringSplit(ss, lines[hasUncSourceWithIndex_s].Data(), " ");

		if(TString(ss[0]).IsDigit() == false) {
			cout<<"   You are trying to assign the following line which not starting from digits to "<< s+1 <<"th systematics sources"<<endl;
			cout<<"\""<<lines[hasUncSourceWithIndex_s]<<"\""<<endl;
			cout<<"   Probable reasons: check kmax or index of systematics sources"<<endl;
			exit(0);
		}

		int indexcorrelation = (TString(ss[0])).Atoi();

		//if(indexcorrelation != (nsyssources-s) ) {
		//	cout<<"Warning: you are reading in "<<nsyssources-s<<"th nuisanse parameter, but this line start with "<<indexcorrelation<<endl;
		//}

		int pdf = 0; 
		if(ss[1]=="lnN") pdf=1;
		else if(ss[1]=="trG") pdf=2;
		else pdf =  (TString(ss[1])).Atoi();

		tmps+= TString::Format("%3s ", ss[1].c_str());

		bool filledThisSource = false;
		for(int p=0; p<ntotprocesses; p++){
			double err = (TString(ss[p+2])).Atof()-1.0;
			tmps+= TString::Format("%7.2f ",err+1.0);
			if(err == 0.) continue;
			//cout<<"delete me: c="<<binnumber[p]-1<<" s="<< subprocess[p]<<endl;
			cms->AddUncertainty(binnumber[p]-1, subprocess[p], err, pdf, indexcorrelation );

			// because when err < 0, AddUncertainty do nothing,  but filledThisSource has been changed to be true
			filledThisSource = true;
		}
		if(!filledThisSource) {
			cout<<"WARNING: The "<< s+1 <<"th source of uncertainties are all 0. "<<endl;
			//FIXME temporarily solution:  when reading a source with all error = 0,  then assign it to be logNormal, error =0,  in UtilsROOT.cc 
			// to incorporate the MinuitFit.  This ad-hoc fix will slow down the toy generation because there are lots of non-use random numebrs and operations ...
			// cms->AddUncertainty(0, 0, 0, 1, indexcorrelation );
			cms->AddUncertainty(0, 0, 0, 1, indexcorrelation );
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

bool ConfigureModel(CountingModel *cms, const char* fileName){
	TString s = ReadFile(fileName);
	return CheckIfDoingShapeAnalysis(cms, s);
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
