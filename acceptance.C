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

void acceptance(string fin="ecut550.lst", bool pinchSeptum=false, double pinch=0.0)
{

  TH1D *hThrow=new TH1D("hThrow","thrown distribution (ev.th); #theta (deg)",200,2,9);
  TH1D *hAccept=new TH1D("hAccept","A(ev.th);#theta (deg)",200,2,9);
  TH1D *hAcceptZ=new TH1D("hAcceptZ","A(th_ztarg); #theta (deg)",200,2,9);
  TH1D *hAcceptNoMSc=new TH1D("hAcceptNoMSc","A(ev.th-MSc); #theta (deg)",200,2,9);
  TH1D *hSigma=new TH1D("hSigma","hist cross section (ev.th); #theta (deg)",200,2,9);
  TH1D *vSigma=new TH1D("vSigma","vtx cross section (ev.th); #theta (deg); average cross section",200,2,9);

  TH1D *zTargIn=new TH1D("zTargIn","events accepted z (Pb) (ev.th);#theta (deg)",200,2,9);
  TH1D *zTargInZ=new TH1D("zTargInZ","events accepted z (Pb) (th_ztarg);#theta (deg)",200,2,9);
  TH1D *zTargOut=new TH1D("zTargOut","events !accepted z (Pb) (ev.th);#theta (deg)",200,2,9);



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

    if(nuclA!=208) continue;
    rate /=nfiles;

    hThrow->Fill(th*rad2deg,rate);
    
    if(th_ztarg==-333){
      zTargOut->Fill(th*rad2deg,rate);
    }else{
      zTargIn->Fill(th*rad2deg,rate);
      zTargInZ->Fill(th_ztarg*rad2deg,rate);
    }
    
    if(xcol==-333) continue;
    if(x_fp_tr==-333) continue;
    if(p_ztarg==-333) continue;
    if(!CollimatorL(xcol, ycol)) continue;
    if( 953.4 - p_ztarg >= 2.2 ) continue;
    
    hAcceptNoMSc->Fill(sampledTh,rate);
    hAccept->Fill(th*rad2deg,rate);
    hAcceptZ->Fill(th_ztarg*rad2deg,rate);

    vSigma->Fill(th*rad2deg,rate*xs);
  }

  //FIXME figure out hSigma

  auto *c1=new TCanvas();
  c1->Divide(3);
  c1->cd(1);
  hThrow->DrawCopy("hist");
  c1->cd(2);
  hAccept->SetLineWidth(3);
  hAccept->DrawCopy("hist");
  hAcceptZ->SetLineColor(2);
  hAcceptZ->DrawCopy("hist && same");
  hAcceptNoMSc->SetLineColor(6);
  hAcceptNoMSc->DrawCopy("hist && same");
  gPad->BuildLegend();
  c1->cd(3);
  TH1D *hd1=(TH1D*)hAccept->Clone("hd1");
  hd1->Divide(hThrow);
  hd1->SetTitle("Acceptance function");
  hd1->DrawCopy("");

  auto *c2=new TCanvas();
  c2->Divide(2);
  c2->cd(1);
  zTargIn->SetLineWidth(3);
  zTargIn->DrawCopy("hist");
  zTargInZ->SetLineColor(2);
  zTargInZ->DrawCopy("hist && same");
  zTargOut->SetLineColor(6);
  zTargOut->DrawCopy("hist && same");
  gPad->BuildLegend();
  c2->cd(2);
  TH1D *hd2=(TH1D*)zTargOut->Clone("hd2");
  hd2->Divide(zTargIn);
  hd2->SetTitle("magenta/blue");
  hd2->DrawCopy("");

  auto *c3=new TCanvas();
  vSigma->Divide(hAccept);
  vSigma->DrawCopy("hist");
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

  tree->SetBranchAddress("ev.Th",    &sampledTh);//this is the sampled theta
  tree->SetBranchAddress("ev.th",    &th);//this is the "lab" theta at vertex
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

  tree->SetBranchAddress("p_ztarg", &p_ztarg);
  tree->SetBranchAddress("ph_ztarg", &ph_ztarg);
  tree->SetBranchAddress("th_ztarg", &th_ztarg);

  tree->SetBranchAddress("p_q1en_tr", &p_q1en_tr);
  tree->SetBranchAddress("p_q1ex_tr", &p_q1ex_tr);
  tree->SetBranchAddress("p_q2ex_tr", &p_q2ex_tr);
  tree->SetBranchAddress("p_dex_tr", &p_dex_tr);
  tree->SetBranchAddress("p_q3ex_tr", &p_q3ex_tr);

  tree->SetBranchAddress("x_fp_tr", &x_fp_tr);

}
