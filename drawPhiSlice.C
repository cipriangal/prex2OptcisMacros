void drawPhiSlice(string fnm){
  TFile *fin=TFile::Open(fnm.c_str(),"READ");
  
  string arm="R";
  if( fin->GetListOfKeys()->Contains("Ldpctr"))
    arm="L";
  else if( !fin->GetListOfKeys()->Contains("Rytg")){
    cout<<"can't figure out arm. quitting"<<endl;
    return;
  }
    
  TTree *t=(TTree*)fin->Get("T");

  auto *c1=new TCanvas();
  c1->Divide(2);
  c1->cd(1);
  t->Draw(Form("%s.gold.th:%s.gold.ph",arm.c_str(),arm.c_str()),Form("abs(%s.gold.ph)<0.03 && abs(%s.gold.th)<0.05",arm.c_str(),arm.c_str()),"colz");

  c1->cd(2);

  //gStyle->SetPalette(kTemperatureMap);
  //gStyle->SetPalette(kVisibleSpectrum);
  int nColors = gStyle->GetNumberOfColors();

  const int nStep=8;
  TH1D *h[nStep];
  float max = 0;
  for(int i=0;i<nStep;i++){
    int histoColor = (float)nColors/nStep * i;

    double low = -0.02 + i*0.04/nStep;
    double hig = -0.02 + (i+1)*0.04/nStep;
    h[i]=new TH1D(Form("h%d",i),Form("%6.3f  < #phi < %6.3f; theta",low,hig),
		  100,-0.05,0.05);
    t->Project(Form("h%d",i),Form("%s.gold.th",arm.c_str()),
	       Form("%s.gold.ph<%f && %s.gold.ph>%f && abs(%s.gold.th)<0.05",
		    arm.c_str(),hig,arm.c_str(),low,arm.c_str()));
    h[i]->SetLineColor(gStyle->GetColorPalette(histoColor));
    h[i]->SetLineWidth(2);
    if(max<h[i]->GetMaximum()) max = h[i]->GetMaximum();
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);

  h[0]->GetYaxis()->SetRangeUser(0,max*1.2);
  h[0]->DrawCopy("");
  for(int i=1;i<nStep;i++)
    h[i]->DrawCopy("SAME");

  gPad->BuildLegend();
}
