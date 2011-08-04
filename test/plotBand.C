{
	/********  usage *********
	 *  
	 root -l PlotUtilities_standalone.cc++    
	 ----it will produce a shared library PlotUtilities_standalone_cc.so
	 root -l plotBand.C
	 *
	 */

	//  use some style template 
	//#include "/data/Projects/FloridaStyle.C"
	//	FloridaStyle();

	gSystem->Load("PlotUtilities_standalone_cc.so");
	gStyle->SetCanvasDefH      (600);
	gStyle->SetCanvasDefW      (800);
	gStyle->SetTitleOffset( 1, "y");

	const int npoints = 5;
	double limits_obs[npoints] = {
		41.74677 ,
		24.58797 ,
		24.73013 ,
		58.37274 ,
		150.18342 , 
	};
	double limits_1fb_obs[npoints] = {
		3.19987 ,
		1.86819 ,
		1.57866 ,
		3.04254 ,
		6.74371 , 
	};
	double bands_1fb[npoints][5] = {
		// -2 sigma, -1 sigma,   mean,     +1 sigma,  +2 sigma
		{ 2.36292  ,  3.10099   , 4.99050   , 6.50240  ,  8.84537},	
		{ 1.34476  ,  1.81795   , 2.92033   , 3.77545  ,  5.19139},
		{ 1.10501  ,  1.33298   , 2.18627   , 2.68424  ,  3.99269},
		{ 2.47239  ,  2.47239   , 3.85914   , 4.53227  ,  6.67285},
		{ 6.21168  ,  6.21168   , 8.17544   , 9.06174  , 12.84551},
	};
	double bands[npoints][5] = {
		// -2 sigma, -1 sigma,   mean,     +1 sigma,  +2 sigma
		{ 43.13062 ,  43.13062 ,  48.80107  , 46.82769 ,  65.03096},
		{ 25.52106 ,  25.52106 ,  28.82003  , 27.57525 ,  38.26959},
		{ 26.90173 ,  26.90173 ,  28.75666  , 26.90173 ,  35.56380},
		{ 65.68063 ,  65.68063 ,  67.36003  , 65.68063 ,  72.34968},
		{170.76523 , 170.76523 , 172.93633  ,170.76523 , 170.76523},
	};
	double xpoints[npoints] = {250, 300, 400, 500, 600};

	bool projectingRLimitLogY = true;
	double projectingXmin = 190, projectingXmax = 610;
	double projectingRLimitYmin = 0.5, projectingRLimitYmax = 200;
	string projectingRLimitXYtitles = ";Higgs mass, m_{H} [GeV/c^{2}]; 95% CL Limit on #sigma/#sigma_{SM} ";
	string ssave= "./limits_35pb";

	double limits_m1s[npoints], limits_m2s[npoints], limits_p1s[npoints], limits_p2s[npoints], limits_mean[npoints], limits_obs[npoints];
	double limits_1fb_m1s[npoints], limits_1fb_m2s[npoints], limits_1fb_p1s[npoints], limits_1fb_p2s[npoints], limits_1fb_mean[npoints], limits_1fb_obs[npoints];
	for(int n=0; n<npoints; n++){
		limits_m2s[n]=bands[n][0];
		limits_m1s[n]=bands[n][1];
		limits_mean[n]=bands[n][2];
		limits_p1s[n]=bands[n][3];
		limits_p2s[n]=bands[n][4];
		limits_1fb_m2s[n]=bands_1fb[n][0];
		limits_1fb_m1s[n]=bands_1fb[n][1];
		limits_1fb_mean[n]=bands_1fb[n][2];
		limits_1fb_p1s[n]=bands_1fb[n][3];
		limits_1fb_p2s[n]=bands_1fb[n][4];
	}

	TPaveText *pt;
	pt = SetTPaveText(0.5, 0.7, 0.8, 0.9);
	pt->AddText("HZZ #rightarrow 2l2#nu");
	PlotWithBelts lb(
			limits_m1s, limits_p1s, limits_m2s, limits_p2s,
			limits_mean, limits_obs, npoints,  
			xpoints, ssave+"_1", pt,
			projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
	lb.plot();
	lb.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: Observed");
	MoveLegend(lb.getLegend(),0.5,0.6);
	lb.getMeanGraph()->SetLineStyle(1);
	lb.getMeanGraph()->SetLineColor(kBlue);
	lb.getMeanGraph()->SetLineWidth(2);
	lb.getObsGraph()->SetLineStyle(2);
	lb.getObsGraph()->SetLineWidth(1);
	lb.getLine()->SetLineColor(kRed);
	PlotWithBelts lb2(
			limits_1fb_m1s, limits_1fb_p1s, limits_1fb_m2s, limits_1fb_p2s,
			limits_1fb_mean, limits_1fb_obs, npoints,  
			xpoints, ssave+"_2", pt,
			projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
	lb2.plot();
	lb2.getObsGraph()->SetMarkerStyle(1);//remove marker
	lb2.drawLegend("95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: mean(nosys)");
	lb2.getMeanGraph()->SetLineStyle(1);
	lb2.getMeanGraph()->SetLineColor(kBlue);
	lb2.getMeanGraph()->SetLineWidth(2);
	lb2.getObsGraph()->SetLineStyle(2);
	lb2.getObsGraph()->SetLineWidth(2);
	lb2.getLine()->SetLineColor(kRed);
	lb.getYellowGraph()->Draw("f");
	lb.getGreenGraph()->Draw("f");
	lb.getMeanGraph()->Draw("l");
	lb.getObsGraph()->SetLineWidth(2);
	lb.getObsGraph()->Draw("l");
	lb2.save();	

	lb.getCanvas()->cd();
	lb2.getYellowGraph()->Draw("f");
	lb2.getGreenGraph()->Draw("f");
	lb2.getMeanGraph()->Draw("l");
	lb2.getObsGraph()->SetLineWidth(2);
	lb2.getObsGraph()->Draw("l");
	lb.save();	
}
