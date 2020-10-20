double pi = TMath::Pi();
double deg2rad = pi/180.;
double rad2deg = 180./pi;

//match to what we used in the analyzer
double amu = 931.494e-3;
double mass_tg = 208*amu;
double mass_e = 0.511e-3;

// event variables
double Q2;
double rate;
double xs;
double thcom, th, ph;
double vx, vy, vz;
int nuclA;
int pid;
double beamp;

double p;
double Asym;
double Am; // scaled by pol
double S;  // sensitivity
double ep;

double p_tg;

// var used for cuts
double xcol, ycol;
double xup1, yup1, xup2, yup2;
double xd1, xd2, xd3, xd4, xd5, xd6, xd7, xd8, xd9;
double yd1, yd2, yd3, yd4, yd5, yd6, yd7, yd8, yd9;

// tar variables
double x_tg, y_tg, z_tg, th_tg, ph_tg;
double th_ztarg_tr, ph_ztarg_tr;

// vdc
double x_vdc_tr, y_vdc_tr, th_vdc_tr, ph_vdc_tr;


//quad,dipole vars
double p_q1en_tr,p_q1ex_tr,p_q2ex_tr,p_dex_tr,p_q3ex_tr;
