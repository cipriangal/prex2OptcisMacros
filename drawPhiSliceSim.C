#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TMath.h>

#include "ryan/CollimatorL.C"
#include "ryan/DownPlane.C"
#include "poleTip.h"

using namespace std;

double E0 = 0.9534;
double angle = 4.8;
double theta0 = angle * deg2rad;

void SetTree(TTree* tree);

void drawPhiSliceSim(string fin="sand.lst", bool pinchSeptum=true, double pinch=0.0)
{

  const int nSlice=8;
  TH2D *tp=new TH2D("tp",";ph_ztarg_tr;th_ztarg_tr",100,-0.03,0.03,100,-0.05,0.05);
  TH1D *p[nSlice];
  float lowLim[nSlice+1];
  for(int i=0;i<nSlice;i++){
    lowLim[i] = -0.025 + i*0.05/nSlice;
    lowLim[i+1] = -0.025 + (i+1)*0.05/nSlice;
    p[i]=new TH1D(Form("p%d",i),Form("%6.3f < #phi < %6.3f",lowLim[i],lowLim[i+1]),100,-0.05,0.05);
  }

  TChain* T = new TChain("T","T");
  SetTree(T);
  
  ifstream ifstr(fin.c_str());
  string fname;
  int nfiles = 0;
  while( ifstr >> fname ){
    T->Add(fname.c_str());
    nfiles++;
  }
  
  int nentries = T->GetEntries();
  for(int ientry=0; ientry<nentries; ientry++){
    T->GetEntry(ientry);

    double costh0 = 1-(Q2/(2*beamp*ep));
    double th_vtx = acos(costh0);
    
    if(CollimatorL(xcol, ycol) && xcol != -333
       && th_ztarg_tr!=-333 && ph_ztarg_tr!=-333){

      double pinchscan = pinch/1000;
      bool scsvdn = DownPlane(xd1,yd1,xd2,yd2,
			      xd3-pinchscan,yd3,xd4-pinchscan,yd4,
			      xd5,yd5,
			      xd6,yd6,xd7,yd7,xd8,yd8,xd9,yd9,1);
      if(!scsvdn && pinchSeptum) continue;

      int found=-2;
      for(int i=0;i<nSlice && found==-2;i++){
	if(ph_ztarg_tr<lowLim[i+1])
	  found=i-1;
      }
      
      if(found>=0 && found<nSlice){
	tp->Fill(ph_ztarg_tr,th_ztarg_tr,rate/5);
	p[found]->Fill(th_ztarg_tr,rate/5);
      }

    }
  }

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2);
  c2->cd(1);
  tp->DrawCopy("colz");
  c2->cd(2);
  int nCls = gStyle->GetNumberOfColors();
  float max=0;
  for(int i=0;i<nSlice;i++){
    int hCl = (float)nCls/nSlice * i;
    
    p[i] -> SetLineWidth(3);
    p[i] -> SetLineColor(gStyle->GetColorPalette(hCl));
    if(max<p[i]->GetMaximum()) max = p[i]->GetMaximum();
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  p[0]->GetYaxis()->SetRangeUser(0,max*1.2);
  p[0]->DrawCopy("hist && c");
  for(int i =1;i<nSlice;i++)
    p[i]->DrawCopy("same&&hist&&c");
  gPad->BuildLegend();
}

void SetTree(TTree* tree)
{
  tree->SetBranchAddress("ev.Q2",    &Q2);
  tree->SetBranchAddress("rate",     &rate);
  tree->SetBranchAddress("ev.thcom", &thcom);

  tree->SetBranchAddress("ev.A",     &Asym);
  tree->SetBranchAddress("ev.Am",    &Am);
  tree->SetBranchAddress("ev.S",     &S);
  tree->SetBranchAddress("ev.beamp", &beamp);
  tree->SetBranchAddress("ev.ep",    &ep);

  //tree->SetBranchAddress("ev.Th",    &th);
  //  tree->SetBranchAddress("ev.th",    &th);
  tree->SetBranchAddress("ev.ph",    &ph);
  tree->SetBranchAddress("ev.pid",   &pid);
  tree->SetBranchAddress("ev.nuclA", &nuclA);
  tree->SetBranchAddress("ev.vx",    &vx);
  tree->SetBranchAddress("ev.vy",    &vy);
  tree->SetBranchAddress("ev.vz",    &vz);
  tree->SetBranchAddress("ev.xs",    &xs);
  tree->SetBranchAddress("ev.p",     &p);

  tree->SetBranchAddress("x_col_tr", &xcol);
  tree->SetBranchAddress("y_col_tr", &ycol);

  tree->SetBranchAddress("x_zup1",   &xup1);
  tree->SetBranchAddress("x_zup2",   &xup2);
  tree->SetBranchAddress("y_zup1",   &yup1);
  tree->SetBranchAddress("y_zup2",   &yup2);

  tree->SetBranchAddress("x_zdown1", &xd1);
  tree->SetBranchAddress("x_zdown2", &xd2);
  tree->SetBranchAddress("x_zdown3", &xd3);
  tree->SetBranchAddress("x_zdown4", &xd4);
  tree->SetBranchAddress("x_zdown5", &xd5);
  tree->SetBranchAddress("x_zdown6", &xd6);
  tree->SetBranchAddress("x_zdown7", &xd7);
  tree->SetBranchAddress("x_zdown8", &xd8);
  tree->SetBranchAddress("x_zdown9", &xd9);

  tree->SetBranchAddress("y_zdown1", &yd1);
  tree->SetBranchAddress("y_zdown2", &yd2);
  tree->SetBranchAddress("y_zdown3", &yd3);
  tree->SetBranchAddress("y_zdown4", &yd4);
  tree->SetBranchAddress("y_zdown5", &yd5);
  tree->SetBranchAddress("y_zdown6", &yd6);
  tree->SetBranchAddress("y_zdown7", &yd7);
  tree->SetBranchAddress("y_zdown8", &yd8);
  tree->SetBranchAddress("y_zdown9", &yd9);

  tree->SetBranchAddress("x_tg",     &x_tg);
  tree->SetBranchAddress("y_tg",     &x_tg);
  tree->SetBranchAddress("z_tg",     &z_tg);
  tree->SetBranchAddress("th_tg",    &th_tg);
  tree->SetBranchAddress("ph_tg",    &ph_tg);

  tree->SetBranchAddress("x_vdc_tr", &x_vdc_tr);
  tree->SetBranchAddress("y_vdc_tr", &y_vdc_tr);
  tree->SetBranchAddress("th_vdc_tr",&th_vdc_tr);
  tree->SetBranchAddress("ph_vdc_tr",&ph_vdc_tr);

  tree->SetBranchAddress("p_ztarg_tr", &p_tg);
  tree->SetBranchAddress("ph_ztarg_tr", &ph_ztarg_tr);
  tree->SetBranchAddress("th_ztarg_tr", &th_ztarg_tr);

  tree->SetBranchAddress("p_q1en_tr", &p_q1en_tr);
  tree->SetBranchAddress("p_q1ex_tr", &p_q1ex_tr);
  tree->SetBranchAddress("p_q2ex_tr", &p_q2ex_tr);
  tree->SetBranchAddress("p_dex_tr", &p_dex_tr);
  tree->SetBranchAddress("p_q3ex_tr", &p_q3ex_tr);

}
