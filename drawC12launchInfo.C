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

void drawC12launchInfo(string fin="sand.lst", bool pinchSeptum=false, double pinch=0.0)
{

  TH1D *hvz[2], *hLth[2], *hLph[2],*hThSep[2],*hPhSep[2],*hQ2[2];
  TH1D *hDotVS[2],*hDotBV[2];
  string hnm[2]={"US","DS"};
  for(int i=0;i<2;i++){
    hvz[i]=new TH1D(Form("hvz_%s",hnm[i].c_str()),Form("%s; vertex z [mm]",hnm[i].c_str())
		    ,120,-0.6,0.6);
    hLth[i]=new TH1D(Form("hLth_%s",hnm[i].c_str()),Form("%s launch; theta [deg]",hnm[i].c_str())
		     ,100,2,8);
    hLph[i]=new TH1D(Form("hLph_%s",hnm[i].c_str()),Form("%s launch; phi [deg]",hnm[i].c_str())
		     ,180,-90,90);
    hThSep[i]=new TH1D(Form("hThSep_%s",hnm[i].c_str()),Form("%s septum entrance; theta [deg]",hnm[i].c_str())
		       ,100,2,8);
    hPhSep[i]=new TH1D(Form("hPhSep_%s",hnm[i].c_str()),Form("%s septum entrance; theta [deg]",hnm[i].c_str())
		       ,180,-90,90);
    hQ2[i]=new TH1D(Form("hQ2_%s",hnm[i].c_str()),Form("%s; Q2 (GeV/c)^2",hnm[i].c_str())
		    ,200,0,0.015);
    hDotVS[i]=new TH1D(Form("hDotVS_%s",hnm[i].c_str()),Form("%s VS;cos(angle)",hnm[i].c_str())
		       ,200,0.995,1.0001);
    hDotBV[i]=new TH1D(Form("hDotBV_%s",hnm[i].c_str()),Form("%s BV;cos(angle)",hnm[i].c_str())
		       ,200,0.995,1.0001);
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
       && th_ztarg_tr!=-333 && ph_ztarg_tr!=-333 && 
       nuclA==12 && p_sen!=-333){

      double pinchscan = pinch/1000;
      bool scsvdn = DownPlane(xd1,yd1,xd2,yd2,
			      xd3-pinchscan,yd3,xd4-pinchscan,yd4,
			      xd5,yd5,
			      xd6,yd6,xd7,yd7,xd8,yd8,xd9,yd9,1);
      if(!scsvdn && pinchSeptum) continue;
      rate /=nfiles;
      
      int ud=0;
      if(vz>0) ud=1;
      hvz[ud]->Fill(vz*1000,rate);
      TVector3 mom(px,py,pz);

      hLth[ud]->Fill(mom.Theta()*rad2deg,rate);
      hLph[ud]->Fill(mom.Phi()*rad2deg,rate);
      hThSep[ud]->Fill(th_sen*rad2deg,rate);
      hPhSep[ud]->Fill(ph_sen*rad2deg,rate);
      hQ2[ud]->Fill(Q2,rate);

      TVector3 pos(vx,vy,8+vz);
      //takes the raster into account
      TVector3 beam=pos.Unit();

      // pre-vertex with MS and raster taken into account
      TVector3 preVbeam(0,0,1);
      preVbeam.SetPhi(b_ph*deg2rad);
      preVbeam.SetTheta(b_th);
      

      TVector3 momSeptum(1,0,0);
      momSeptum.SetMag(p_sen);
      momSeptum.SetPhi(ph_sen);
      momSeptum.SetTheta(th_sen);
      hDotVS[ud]->Fill(mom*momSeptum/(mom.Mag()*momSeptum.Mag()),rate);
      hDotBV[ud]->Fill(beam*preVbeam/(beam.Mag()*preVbeam.Mag()),rate);
    }
  }

  hvz[1]->SetLineColor(2);
  hLth[1]->SetLineColor(2);
  hThSep[1]->SetLineColor(2);
  hPhSep[1]->SetLineColor(2);
  hLph[1]->SetLineColor(2);
  hQ2[1]->SetLineColor(2);
  hDotVS[1]->SetLineColor(2);
  hDotBV[1]->SetLineColor(2);

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->Divide(2,2);
  c2->cd(1);
  hvz[0]->DrawCopy("hist");
  hvz[1]->DrawCopy("hist && samej");
  gPad->BuildLegend();
  c2->cd(2);
  hQ2[0]->DrawCopy("hist");
  hQ2[1]->DrawCopy("hist && samej");
  gPad->BuildLegend();
  c2->cd(3);
  hDotVS[1]->DrawCopy("hist");
  hDotVS[0]->DrawCopy("hist && same");
  gPad->SetLogy(1);
  gPad->BuildLegend();
  c2->cd(4);
  hDotBV[1]->DrawCopy("hist");
  hDotBV[0]->DrawCopy("hist && same");
  gPad->SetLogy(1);
  gPad->BuildLegend();

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->Divide(2,2);
  c1->cd(1);
  hLth[1]->DrawCopy("hist");
  hLth[0]->DrawCopy("hist && samej");
  gPad->BuildLegend();
  c1->cd(2);
  hThSep[0]->DrawCopy("hist");
  hThSep[1]->DrawCopy("hist && samej");
  gPad->BuildLegend();
  c1->cd(3);
  hLph[0]->DrawCopy("hist");
  hLph[1]->DrawCopy("hist && samej");
  gPad->BuildLegend();
  c1->cd(4);
  hPhSep[0]->DrawCopy("hist");
  hPhSep[1]->DrawCopy("hist && samej");
  gPad->BuildLegend();

  TCanvas *c3=new TCanvas("c3","c3");
  c3->Divide(2);
  c3->cd(1);
  hDotBV[1]->DrawCopy("hist");
  hDotVS[0]->DrawCopy("hist&&same");
  gPad->SetLogy(1);
  gPad->BuildLegend();
  c3->cd(2);
  hDotBV[0]->DrawCopy("hist");
  hDotVS[1]->DrawCopy("hist&&same");
  gPad->SetLogy(1);
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
  tree->SetBranchAddress("ev.px",     &px);
  tree->SetBranchAddress("ev.py",     &py);
  tree->SetBranchAddress("ev.pz",     &pz);

  tree->SetBranchAddress("bm.th",&b_th);
  tree->SetBranchAddress("bm.ph",&b_ph);

  tree->SetBranchAddress("th_sen", &th_sen);
  tree->SetBranchAddress("ph_sen", &ph_sen);
  tree->SetBranchAddress("p_sen", &p_sen);

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
