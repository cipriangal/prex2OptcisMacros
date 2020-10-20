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
#include "poleTip.h"

using namespace std;

double E0 = 0.9534;
double angle = 4.8;
double theta0 = angle * deg2rad;

void SetTree(TTree* tree);

void poleTip(string fin="nom.lst")
{

  string hnm[5]={"q1in","q1out","q2out","dout","q3out"};
  TH1D *p[5];
  for(int i=0;i<5;i++)
    p[i]=new TH1D(Form("p%s",hnm[i].c_str()),Form("momentum %s",hnm[i].c_str()),100,860,955);

  TChain* T = new TChain("T","T");
  SetTree(T);
  
  ifstream ifstr(fin.c_str());
  string fname;
  int nfiles = 0;
  while( ifstr >> fname )
    {
      T->Add(fname.c_str());
      nfiles++;
    }

  int nentries = T->GetEntries();
  for(int ientry=0; ientry<nentries; ientry++){
    T->GetEntry(ientry);
    double costh0 = 1-(Q2/(2*beamp*ep));
    double th_vtx = acos(costh0);
    
    if(CollimatorL(xcol, ycol) && xcol != -333){

      if(p_q1en_tr!=-333)
	p[0]->Fill(p_q1en_tr,rate);

      if(p_q1ex_tr!=-333)
	p[1]->Fill(p_q1ex_tr,rate);

      if(p_q2ex_tr!=-333)
	p[2]->Fill(p_q2ex_tr,rate);

      if(p_dex_tr!=-333)
	p[3]->Fill(p_dex_tr,rate);

      if(p_q3ex_tr!=-333)
	p[4]->Fill(p_q3ex_tr,rate);
    }
  }

  TCanvas *c2 = new TCanvas("c2", "c2");
  p[0]->SetLineColor(1);
  cout<<"Integra;ls"<<endl;
  p[0]->DrawCopy();
  cout<<p[0]->Integral()/10<<endl;
  for(int i=1;i<5;i++){
    p[i]->SetLineColor(i+1);
    p[i]->DrawCopy("same");
    cout<<"\t"<<hnm[i]<<" "<<p[i]->Integral()/10<<endl;
  }

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
  tree->SetBranchAddress("ph_ztarg_tr", &thisph);
  tree->SetBranchAddress("th_ztarg_tr", &thisth);

  tree->SetBranchAddress("p_q1en_tr", &p_q1en_tr);
  tree->SetBranchAddress("p_q1ex_tr", &p_q1ex_tr);
  tree->SetBranchAddress("p_q2ex_tr", &p_q2ex_tr);
  tree->SetBranchAddress("p_dex_tr", &p_dex_tr);
  tree->SetBranchAddress("p_q3ex_tr", &p_q3ex_tr);

}
