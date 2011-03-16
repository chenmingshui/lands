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
				if(lines[k].BeginsWith("Observation") or lines[k].BeginsWith("Observation")){
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
					if(TString(ss1[0]).Atoi()!=u) continue; // FIXME  need to adapt to names 
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
			break;

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
		if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==true){
			cout<<"WARNING:  there are two lines beginning with \"bin\" in your card. We will take the first one"<<endl;
		}
		if(lines[l].BeginsWith("bin ") && hasLineStartWithBin==false){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");

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
		if(lines[l].BeginsWith("bins ") or lines[l].BeginsWith("binname") ){
			vector<string>	ss; 
			ss.clear();
			StringSplit(ss, lines[l].Data(), " ");
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
	nsyssources = 0;
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
		else pdf =  (TString(ss[1])).Atoi();

		tmps+= TString::Format("%3s ", ss[1].c_str());
		if(pdf==typeGamma && ss[1]!="gmM") tmps+= TString::Format("%8s ", ss[2].c_str());
		else		   tmps+= "         ";

		bool filledThisSource = false;
		for(int p=0; p<ntotprocesses; p++){
			double err, rho;
			if(pdf==typeLogNormal){
				if(ss[p+2]=="-") {
					tmps+= "    -   ";
					continue;
				}
			       	err= (TString(ss[p+2])).Atof()-1.0;
				tmps+= TString::Format("%7.2f ",err+1.0);
				if(err < -1) {
					cout<<"Kappa in lognormal  can't be negative, please check your input card at "<<s+1<<"th source, "<<p+2<<"th entry"<<endl;
					exit(0);
				}
				if(err == 0.) continue;
			}
			else if(pdf==typeTruncatedGaussian){
				if(ss[p+2]=="-"){
					tmps+= "    -   ";
				       	continue;
				}
			       	err= (TString(ss[p+2])).Atof();
				tmps+= TString::Format("%7.2f ",err);
				if(err == 0.) continue;
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
			if(pdf==typeLogNormal||pdf==typeTruncatedGaussian)cms->AddUncertainty(binnumber[p]-1, subprocess[p], err, pdf, indexcorrelation );
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
