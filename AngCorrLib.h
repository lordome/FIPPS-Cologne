#include "Riostream.h"
#include "TRint.h"
#include "TGClient.h"
#include "TString.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TBox.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "Math/IFunction.h"
#include "TLatex.h"
#include "TFrame.h"
#include <cmath>
#include "Math/SpecFuncMathMore.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TMarker.h"
#include "TLine.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TGraphAsymmErrors.h"
#include "TStatistic.h"
#include "TRandom3.h"

#include "FIPPS_Softs.h"

using namespace  std;

bool fUseMonteCarlo=false;
Int_t fNMonteCarlo=10000;

TCanvas *fCanvas = nullptr;

array<TGraphErrors *,3> gAngCorr;
array<TGraphErrors *,3> gAngCorrNorm;
array<TGraphErrors *,3> gAngCorrCorr;

TGraph *gDummy=nullptr;
TGraph *gA2A4=nullptr;
TGraphAsymmErrors *fCurrentPoint = nullptr;
TGraphErrors *fCurrentTheoPoint = nullptr;
double BestDelta;

array < map < Int_t, TH2*>, 3 > vListOfMatrices;
array < map < Int_t, TH2*>, 3 > vListOfNormMatrices;

pair<Float_t, Float_t> GatePx;
pair<Float_t, Float_t> GatePy;

pair<Float_t, Float_t> GetPx() {return GatePx;}
pair<Float_t, Float_t> GetPy() {return GatePy;}

pair<Float_t, Float_t> GateBx;
pair<Float_t, Float_t> GateBy;

pair<Float_t, Float_t> GetBx() {return GateBx;}
pair<Float_t, Float_t> GetBy() {return GateBy;}

pair<Float_t, Float_t> GateB2x;
pair<Float_t, Float_t> GateB2y;

pair<Float_t, Float_t> GetB2x() {return GateB2x;}
pair<Float_t, Float_t> GetB2y() {return GateB2y;}

vector<double> vminima;

array<TF1*,3> fAngCorFuncFit;
array<TGraphAsymmErrors*,3> fAngCorFuncFit_Err;

array<TF1*,3> fAngCorFuncFitTheo;

ROOT::Math::Minimizer *fMinimum = nullptr;

TF1 *fAngCorrFunc = nullptr;
TGraphAsymmErrors *fAngCorrFuncErr = nullptr;

TF1 *fTheoAngCorrFunc = nullptr;

TF1 *fChi2Func = nullptr;
TF1 *fChi2FuncDeriv = nullptr;
TF1 *fChi2FuncDeriv2 = nullptr;

TH2F *fChi2Func2D = nullptr;


TBox *boxPxPy  = nullptr;
TBox *boxPxBy  = nullptr;
TBox *boxBxPy  = nullptr;
TBox *boxBxBy  = nullptr;
TBox *boxBxB2y = nullptr;
TBox *boxB2xBy = nullptr;
TBox *boxB2xB2y = nullptr;

TString fCurrentMode;
//float Q2 = 1.;
//float Q4 = 1.;

enum DetType { kMode_FF, kMode_FI, kMode_II};
TString DetTypesString[3] = {"FIPPS-FIPPS","FIPPS-IFIN","IFIN-IFIN"};

// pour chaque type de combinaison de détecteur, on a Q2 +- EQ2 et Q4 +- EQ4
array <pair<pair<double,double>,pair<double,double>>,3> fQValues;

int ftwoI1=0;
int ftwoI2=0;
int ftwoI3=0;
double fmix1=0.;
double fmix2=0.;

int ftwoL1=0;
int ftwoL1p=0;
int ftwoL2=0;
int ftwoL2p=0;
int kmax=0;
std::vector<double> Akks;

double fA0Exp;
double fA0Exp_errlow;
double fA0Exp_errhigh;
double fA2Exp;
double fA2Exp_errlow;
double fA2Exp_errhigh;
double fA4Exp;
double fA4Exp_errlow;
double fA4Exp_errhigh;

// Angle, Mode, N
map <Int_t , map<Int_t,Int_t > > fMapOfAngles;
map <Int_t , map<Int_t,Int_t > > fMapOfAnglesNorm;

Int_t fFitMixingMode=1;

void SetMonteCarlo(bool on=true, Int_t N = 10000) {fUseMonteCarlo = on; fNMonteCarlo=N;}
void ShowAngles();
void Fit(double xmin=0., double xmax=180., bool same=false);
Double_t AngCorrFuncion(Double_t*xx,Double_t*pp);
void SetTransition(int twoI1, int twoI2, int twoI3, double mix1=0., double mix2=0.);
double Fk(int twoL1, int twoL1p, int twoIi, int twoI, int k);
void Eval_Ak();
double Evaluate(double *x,double *parameters);
void UpdateA2A4();
void UpdateA2A42D();
double EvalChi2(double *x,double *parameters);
double EvalChi2_2D(double atandelta1,double atandelta2);
double EvalChi2Deriv(double *x,double */*parameters*/);
double EvalChi2Deriv2(double *x,double */*parameters*/);

void PrintQValues() {
    cout << "Q values:"<<endl;
    for(int i=0 ; i<3 ; i++) {
        cout<<left<<setw(12)<<DetTypesString[i]<<": Q2 = " << setw(10) << fQValues.at(i).first.first << " +- " << setw(10) << fQValues.at(i).first.second<< endl;
        cout<<left<<setw(12)<<""               <<": Q4 = " << setw(10) << fQValues.at(i).second.first << " +- " << setw(10) << fQValues.at(i).second.second<< endl;
    }
}

void SetQ2Q4(double q2, double q4, double eq2, double eq4, int mode) {
    if(mode<0 || mode>2) {
        ERR_MESS << "mode should be 0 for FIPPS-FIPPS, 1 for FIPPS-IFIN, and 2 for IFIN-IFIN" << ENDL;
    }
    fQValues.at(mode).first.first = q2;
    fQValues.at(mode).first.second = eq2;
    fQValues.at(mode).second.first = q4;
    fQValues.at(mode).second.second = eq4;
    PrintQValues();
}

void ReadQValues() {
    // Read Q values
    ifstream qval_file("qvalues.txt");
    // if not exists, create a dummy one
    if(!qval_file) {
        ofstream qval_file_o("qvalues.txt");
        qval_file_o << left << setw(12) << "Det types  : " << setw(10) << "A2" << setw(10) << "A2 unc" << setw(10) << "A4" << setw(10) << "A4 unc" << endl;
        qval_file_o << left << setw(12) << "FIPPS-FIPPS: " << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << endl;
        qval_file_o << left << setw(12) << "FIPPS-IFIN : " << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << endl;
        qval_file_o << left << setw(12) << "IFIN-IFIN  : " << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << setw(10) << 1. << endl;
        qval_file_o.close();
        qval_file.open("qvalues.txt");
    }
    if(qval_file) {
        string line;
        getline(qval_file,line);
        for(int i=0 ; i<3 ; i++) {
            getline(qval_file,line);
            TString Line = line;
            TObjArray *arr = Line.Tokenize(" ");
            fQValues.at(i).first.first   = ((TString)arr->At(arr->GetEntries()-4)->GetName()).Atof();
            fQValues.at(i).first.second  = ((TString)arr->At(arr->GetEntries()-3)->GetName()).Atof();
            fQValues.at(i).second.first  = ((TString)arr->At(arr->GetEntries()-2)->GetName()).Atof();
            fQValues.at(i).second.second = ((TString)arr->At(arr->GetEntries()-1)->GetName()).Atof();
            delete arr;
        }
        qval_file.close();
    }

    cout<<"Q values read from: qvalues.txt"<<endl;
    PrintQValues();
}

void LoadFile(TString FileName, TString Cond = "NoCond", bool Symetrize=false) {
    if(!gSystem->IsFileInIncludePath(FileName)) {
        ERR_MESS << FileName << " not found " << ENDL;
        return;
    }

    ReadQValues();

    TFile *file = TFile::Open(FileName);
    TDirectoryFile *dir1 = (TDirectoryFile*) file->Get("RelAngles");
    if(dir1 == nullptr) {
        ERR_MESS << FileName << " does not contain angular correlation matrices" << ENDL;
        return;
    }

    TDirectoryFile *dir2 = (TDirectoryFile*) dir1->Get(Cond);
    if(dir2 == nullptr) {
        ERR_MESS << "No matrices for Gate: " << Cond << ENDL;
        dir1->ls();
        return;
    }

    for(int i=0 ; i<3 ; i++) {
        vListOfMatrices[i].clear();
        vListOfNormMatrices[i].clear();
    }


    // Loop on the modes
    for(int i=0 ; i<dir2->GetListOfKeys()->GetEntries() ; i++) {
        TString Name = dir2->GetListOfKeys()->At(i)->GetName();
        Int_t Mode;
        if(Name==DetTypesString[0]) Mode=0;
        else if(Name==DetTypesString[1]) Mode=1;
        else if(Name==DetTypesString[2]) Mode=2;
        else {
            cout << "OUPS, matrices mode: " << Name << " unkown, ignored" << endl;
            continue;
        }

        TDirectoryFile *dir3 = (TDirectoryFile*) dir2->Get(Name);
        if(dir3 == nullptr) {
            ERR_MESS << "No matrices for Gate: " << Cond << " and mode" << Name << ENDL;
            dir2->ls();
            continue;
        }

        // Loop on the matrices
        for(int j=0 ; j<dir3->GetListOfKeys()->GetEntries() ; j++) {
            TString Name2 = dir3->GetListOfKeys()->At(j)->GetName();

            TObject *o = dir3->Get(Name2);

            if(!o->InheritsFrom(TH2::Class())) continue;
            if(((TH2*)o)->GetEntries()==0) continue;

            TObjArray *arr = Name2.Tokenize("_");

            int Angle = ((TString)arr->At(3)->GetName()).Atoi();
            if(Angle <= 28) continue;

            if(Symetrize) {
                if(Mode==0) {
                    if(Angle==180) Angle=0;
                    if(Angle==135) Angle=45;
                }
                if(Mode==1) {
                    if(Angle==135) Angle=45;
                    if(Angle==120) Angle=60;
                }
                if(Mode==2) {
                    if(Angle==180) Angle=0;
                    if(Angle==120) Angle=60;
                }
            }

            bool Norm = Name2.Contains("Norm");

            if(Norm) {
                if(vListOfNormMatrices[Mode].count(Angle)) vListOfNormMatrices[Mode][Angle]->Add((TH2*)o);
                else vListOfNormMatrices[Mode][Angle] = (TH2*)o;
            }
            else {
                if(vListOfMatrices[Mode].count(Angle)) vListOfMatrices[Mode][Angle]->Add((TH2*)o);
                else vListOfMatrices[Mode][Angle] = (TH2*)o;
            }
        }
    }

    ShowAngles();
}

void ShowAngles() {

    for(int imode=0 ; imode<3 ; imode++) {
        if(vListOfMatrices[imode].size()) {
            INFO_MESS << DetTypesString[imode] << " Angles: " << vListOfMatrices[imode].size() << ENDL;
            INFO_MESS << "Available angles: ";
            for(auto i : vListOfMatrices[imode]) cout<< setw(5)<<i.first;
            cout<<ENDL;
            INFO_MESS << DetTypesString[imode] << " Normalisation Angles: " << vListOfNormMatrices[imode].size() << ENDL;
            INFO_MESS << "Available angles: ";
            for(auto i : vListOfNormMatrices[imode]) cout<< setw(5)<<i.first;
            cout<<ENDL;
        }
    }
}


void UpdateGates(bool autorange = false) {
    if(GatePx.first ==0.) {
        ERR_MESS << "You need first to define the gates using: SetGate(Px,w_Px,Py,w_Py)" << ENDL;
        return;
    }

    if(fCanvas==nullptr) {
        ERR_MESS << "No matrix loaded, use PlotAngle(Theta)" << ENDL;
        return;
    }

    fCanvas->cd(1);
    TH2 *hist = nullptr;
    for(int i=0 ; i<gPad->GetListOfPrimitives()->GetEntries() ; i++) {
        TObject *o = gPad->GetListOfPrimitives()->At(i);
        if(o->InheritsFrom(TH2::Class())) hist = ((TH2*)o);
        if(o->InheritsFrom(TBox::Class())) {
            gPad->GetListOfPrimitives()->Remove(o);
            i--;
        }
    }

    if(hist == nullptr) {
        ERR_MESS << "No matrix loaded, use PlotAngle(Theta)" << ENDL;
        return;
    }

    if(autorange) {
        hist->GetXaxis()->SetRangeUser(GatePx.first-5*GatePx.second,GatePx.first+5*GatePx.second);
        hist->GetYaxis()->SetRangeUser(GatePy.first-5*GatePy.second,GatePy.first+5*GatePy.second);
    }

    delete boxPxPy;boxPxPy = nullptr;
    delete boxPxBy;boxPxBy = nullptr;
    delete boxBxPy;boxBxPy = nullptr;
    delete boxBxBy;boxBxBy = nullptr;
    delete boxBxB2y;boxBxB2y = nullptr;
    delete boxB2xBy;boxB2xBy = nullptr;
    delete boxB2xB2y;boxB2xB2y = nullptr;

    boxPxPy = new TBox(GatePx.first-GatePx.second*0.5,GatePy.first-GatePy.second*0.5,GatePx.first+GatePx.second*0.5,GatePy.first+GatePy.second*0.5);
    boxPxPy->SetLineColor(kBlack);
    boxPxPy->SetFillStyle(3001);
    boxPxPy->SetFillColorAlpha(kBlack,0.75);
    boxPxPy->SetLineWidth(3);
    boxPxPy->Draw("same");

    if(GateBx.first > 0.) {
        boxPxBy = new TBox(GatePx.first-GatePx.second*0.5,GateBy.first-GateBy.second*0.5,GatePx.first+GatePx.second*0.5,GateBy.first+GateBy.second*0.5);
        boxPxBy->SetLineColor(kRed);
        boxPxBy->SetFillStyle(3001);
        boxPxBy->SetFillColorAlpha(kRed,0.75);
        boxPxBy->SetLineWidth(3);
        boxPxBy->Draw("same");

        boxBxPy = new TBox(GateBx.first-GateBx.second*0.5,GatePy.first-GatePy.second*0.5,GateBx.first+GateBx.second*0.5,GatePy.first+GatePy.second*0.5);
        boxBxPy->SetLineColor(kRed);
        boxBxPy->SetFillStyle(3001);
        boxBxPy->SetFillColorAlpha(kRed,0.75);
        boxBxPy->SetLineWidth(3);
        boxBxPy->Draw("same");

        boxBxBy = new TBox(GateBx.first-GateBx.second*0.5,GateBy.first-GateBy.second*0.5,GateBx.first+GateBx.second*0.5,GateBy.first+GateBy.second*0.5);
        boxBxBy->SetLineColor(kCyan);
        boxBxBy->SetFillStyle(3001);
        boxBxBy->SetFillColorAlpha(kCyan,0.75);
        boxBxBy->SetLineWidth(3);
        boxBxBy->Draw("same");
    }
    if(GateB2x.first > 0.) {
        boxBxB2y = new TBox(GateBx.first-GateBx.second*0.5,GateB2y.first-GateB2y.second*0.5,GateBx.first+GateBx.second*0.5,GateB2y.first+GateB2y.second*0.5);
        boxBxB2y->SetLineColor(kRed);
        boxBxB2y->SetFillStyle(3001);
        boxBxB2y->SetFillColorAlpha(kRed,0.75);
        boxBxB2y->SetLineWidth(3);
        boxBxB2y->Draw("same");

        boxB2xBy = new TBox(GateB2x.first-GateB2x.second*0.5,GateBy.first-GateBy.second*0.5,GateB2x.first+GateB2x.second*0.5,GateBy.first+GateBy.second*0.5);
        boxB2xBy->SetLineColor(kRed);
        boxB2xBy->SetFillStyle(3001);
        boxB2xBy->SetFillColorAlpha(kRed,0.75);
        boxB2xBy->SetLineWidth(3);
        boxB2xBy->Draw("same");

        boxB2xB2y = new TBox(GateB2x.first-GateB2x.second*0.5,GateB2y.first-GateB2y.second*0.5,GateB2x.first+GateB2x.second*0.5,GateB2y.first+GateB2y.second*0.5);
        boxB2xB2y->SetLineColor(kCyan);
        boxB2xB2y->SetFillStyle(3001);
        boxB2xB2y->SetFillColorAlpha(kCyan,0.75);
        boxB2xB2y->SetLineWidth(3);
        boxB2xB2y->Draw("same");
    }
}

void InitCanvas(TString Mode)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    fCanvas = new TCanvas("Angular_Correlations","Angular Correlations",1600,1000);
    fCanvas->Divide(2,2,0.001,0.001);
    fCanvas->GetPad(1)->SetMargin(0.128677,0.0218672,0.101717,0.0302172 );
    fCanvas->GetPad(2)->SetMargin(0.128677,0.0218672,0.101717,0.0302172 );
    fCanvas->GetPad(3)->SetMargin(0.128677,0.0218672,0.101717,0.0302172 );
    fCanvas->GetPad(4)->SetMargin(0.128677,0.0218672,0.101717,0.0302172 );

    fCanvas->GetPad(3)->Divide(1,2,0.001,0.001);
    if(Mode.CountChar('1')==3){
        fCanvas->GetPad(3)->cd(1)->Divide(3,1,0.001,0.001);
        for(int i=0 ; i<3 ; i++) {
            fCanvas->GetPad(3)->cd(1)->cd(i+1)->SetMargin(0.14369,0.0316819,0.101619,0.00744921);
        }
        fCanvas->GetPad(3)->cd(2)->SetMargin(0.0448321,0.0120952,0.101485,0.096212);
    }
    else if(Mode.CountChar('1')==2) {
        fCanvas->GetPad(3)->cd(1)->Divide(2,1,0.001,0.001);
        for(int i=0 ; i<2 ; i++) {
            fCanvas->GetPad(3)->cd(1)->cd(i+1)->SetMargin(0.14369,0.0316819,0.101619,0.00744921);
        }
        fCanvas->GetPad(3)->cd(2)->SetMargin(0.0448321,0.0120952,0.101485,0.096212);
    }
    else {
        fCanvas->GetPad(3)->cd(1)->SetMargin(0.14369,0.0316819,0.101619,0.00744921);
        fCanvas->GetPad(3)->cd(2)->SetMargin(0.0448321,0.0120952,0.101485,0.096212);
    }

    fCanvas->GetPad(1)->SetLogz();

    fCanvas->ToggleEventStatus();
}

void PlotAllAngle(TString Mode="111", bool Norm = false)
{
    fCurrentMode = Mode;

    TH2 *Fullhist = nullptr;

    if(Mode.Length()!=3) {
        ERR_MESS << "Mode needs to be a 3 digits string (111,011...) corresponding to FIPPS-FIPPS, FIPPS-IFIN, IFIN-IFIN" << ENDL;
        return;
    }

    for(int imode=0 ; imode<3 ; imode++) {
        if(imode==0 && (Mode[0] != '1')) continue;
        if(imode==1 && (Mode[1] != '1')) continue;
        if(imode==2 && (Mode[2] != '1')) continue;

        if(!Norm) {
            for(auto &i : vListOfMatrices[imode]) {
                TH2 *hist = i.second;
                if(!Fullhist) Fullhist = (TH2*) hist->Clone("All_anles");
                else Fullhist->Add(hist);
            }
        }
        else {
            for(auto &i : vListOfNormMatrices[imode]) {
                TH2 *hist = i.second;
                if(!Fullhist) Fullhist = (TH2*) hist->Clone("All_anles");
                else Fullhist->Add(hist);
            }
        }
    }

    if(Fullhist == nullptr) {
        ERR_MESS << "No matrices found for mode " << Mode << endl;
        return;
    }
    delete fCanvas;
    InitCanvas(Mode);

    fCanvas->cd(1);
    Fullhist->Draw("col");
    Fullhist->GetXaxis()->SetTitle("E#gamma 1 (keV)");
    Fullhist->GetYaxis()->SetTitle("E#gamma 2 (keV)");
    Fullhist->GetXaxis()->SetTitleSize(0.05);
    Fullhist->GetXaxis()->SetLabelSize(0.05);
    Fullhist->GetYaxis()->SetTitleSize(0.05);
    Fullhist->GetYaxis()->SetLabelSize(0.05);
    Fullhist->GetXaxis()->CenterTitle();
    Fullhist->GetYaxis()->CenterTitle();


    if(GatePx.first > 0.) {
        UpdateGates(true);
    }
}

void PlotAngle(Int_t Angle, TString Mode="111", bool Norm = false)
{
    fCurrentMode = Mode;

    TH2 *Fullhist = nullptr;

    if(Mode.Length()!=3) {
        ERR_MESS << "Mode needs to be a 3 digits string (111,011...) corresponding to FIPPS-FIPPS, FIPPS-IFIN, IFIN-IFIN" << ENDL;
        return;
    }

    for(int imode=0 ; imode<3 ; imode++) {
        if(imode==0 && (Mode[0] != '1')) continue;
        if(imode==1 && (Mode[1] != '1')) continue;
        if(imode==2 && (Mode[2] != '1')) continue;

        if(!Norm) {
            for(auto &i : vListOfMatrices[imode]) {
                if(i.first != Angle) continue;
                TH2 *hist = i.second;
                if(!Fullhist) Fullhist = (TH2*) hist->Clone(Form("Angle_%d",Angle));
                else Fullhist->Add(hist);
            }
        }
        else {
            for(auto &i : vListOfNormMatrices[imode]) {
                if(i.first !=Angle) continue;
                TH2 *hist = i.second;
                if(!Fullhist) Fullhist = (TH2*) hist->Clone(Form("Angle_%d_Norm",Angle));
                else Fullhist->Add(hist);
            }
        }
    }

    if(Fullhist==nullptr) {
        if(Norm) ERR_MESS << "Normalized matrix for angle Angle " << Angle << " is not present in this dataset" << endl;
        else ERR_MESS << "matrix for angle Angle " << Angle << " is not present in this dataset" << endl;
        return;
    }

    delete fCanvas;
    InitCanvas(Mode);

    fCanvas->cd(1);

    Fullhist->Draw("col");
    Fullhist->GetXaxis()->SetTitle("E#gamma 1 (keV)");
    Fullhist->GetYaxis()->SetTitle("E#gamma 2 (keV)");
    Fullhist->GetXaxis()->SetTitleSize(0.05);
    Fullhist->GetXaxis()->SetLabelSize(0.05);
    Fullhist->GetYaxis()->SetTitleSize(0.05);
    Fullhist->GetYaxis()->SetLabelSize(0.05);
    Fullhist->GetXaxis()->CenterTitle();
    Fullhist->GetYaxis()->CenterTitle();


    if(GatePx.first > 0.) {
        UpdateGates(true);
    }
}


void SetGate(Float_t Px, Float_t w_Px, Float_t Py, Float_t w_Py, bool inverse = false) {
    GatePx.first = Px;
    GatePx.second = w_Px;

    GatePy.first = Py;
    GatePy.second = w_Py;

    GateBx.first = GatePx.first-1.5*GatePx.second;
    if(inverse) GateBx.first = GatePx.first+1.5*GatePx.second;
    GateBx.second = GatePx.second;

    GateBy.first = GatePy.first-1.5*GatePy.second;
    if(inverse) GateBy.first = GatePy.first+1.5*GatePy.second;
    GateBy.second = GatePy.second;

    GateB2x.first = GatePx.first+GatePy.first-GateBy.first;
    GateB2x.second = GatePx.second = GatePx.second;

    GateB2y.first = GatePx.first+GatePy.first-GateBx.first;
    GateB2y.second = GatePy.second;

    UpdateGates(true);

    // Check and adapt if necessary the widths of the gates to correspond to the closests bins

    vector<TBox *> boxes{boxPxPy,boxPxBy,boxBxPy,boxBxBy,boxBxB2y,boxB2xBy,boxB2xB2y};
    TH1 *hist_test = vListOfMatrices.front().begin()->second;
    for(auto box: boxes) {
        box->SetX1(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetXaxis()->FindBin(box->GetX1())));
        box->SetX2(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetXaxis()->FindBin(box->GetX2())));
        box->SetY1(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetYaxis()->FindBin(box->GetY1())));
        box->SetY2(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetYaxis()->FindBin(box->GetY2())));
    }

    fCanvas->GetPad(1)->Modified();
    fCanvas->GetPad(1)->Update();
}

void SetGateX(Float_t Px, Float_t w_Px, bool inverse = false) {
    GatePx.first = Px;
    GatePx.second = w_Px;

    GateBx.first = GatePx.first-1.5*GatePx.second;
    if(inverse) GateBx.first = GatePx.first+1.5*GatePx.second;
    GateBx.second = GatePx.second;

    GateB2x.first = GatePx.first+GatePy.first-GateBy.first;
    GateB2x.second = GatePx.second = GatePx.second;

    GateB2y.first = GatePx.first+GatePy.first-GateBx.first;
    GateB2y.second = GatePy.second;

    UpdateGates();
}

void SetGateY(Float_t Py, Float_t w_Py, bool inverse = false) {

    GatePy.first = Py;
    GatePy.second = w_Py;

    GateBy.first = GatePy.first-1.5*GatePy.second;
    if(inverse) GateBy.first = GatePy.first+1.5*GatePy.second;
    GateBy.second = GatePy.second;

    GateB2x.first = GatePx.first+GatePy.first-GateBy.first;
    GateB2x.second = GatePx.second = GatePx.second;

    GateB2y.first = GatePx.first+GatePy.first-GateBx.first;
    GateB2y.second = GatePy.second;

    UpdateGates();
}


void SetBackGround(Float_t Bx, Float_t By) {
    GateBx.first = Bx;
    GateBx.second = GatePx.second;

    GateBy.first = By;
    GateBy.second = GatePy.second;

    GateB2x.first = GatePx.first+GatePy.first-GateBy.first;
    GateB2x.second = GatePx.second = GatePx.second;

    GateB2y.first = GatePx.first+GatePy.first-GateBx.first;
    GateB2y.second = GatePy.second;

    UpdateGates();
}

void SetBackGroundX(Float_t Bx) {
    GateBx.first = Bx;
    GateBx.second = GatePx.second;

    GateB2y.first = GatePx.first+GatePy.first-GateBx.first;
    GateB2y.second = GatePy.second;

    UpdateGates();
}

void SetBackGroundY(Float_t By) {
    GateBy.first = By;
    GateBy.second = GatePy.second;

    GateB2x.first = GatePx.first+GatePy.first-GateBy.first;
    GateB2x.second = GatePx.second = GatePx.second;

    UpdateGates();
}

void SetAdvancedBackGround(Float_t B2x, Float_t B2y) {
    GateB2x.first = B2x;
    GateB2x.second = GatePx.second;

    GateB2y.first = B2y;
    GateB2y.second = GatePy.second;

    UpdateGates();
}

void SetAdvancedBackGroundX(Float_t B2x) {
    GateB2x.first = B2x;
    GateB2x.second = GatePx.second;

    UpdateGates();
}

void SetAdvancedBackGroundY(Float_t B2y) {
    GateB2y.first = B2y;
    GateB2y.second = GatePy.second;

    UpdateGates();
}

bool ProcessProject(Int_t DoNorm, bool DoStdBgd, bool DoAdvBgd, bool NormTo90=true, bool Print=true) {

    if(!vListOfMatrices.size()) {
        ERR_MESS << "Correlation angles matrices not loaded, use LoadFile(FileName,Cond)" << ENDL;
        return false;
    }
    if(DoNorm && !vListOfNormMatrices.size()) {
        ERR_MESS << "No normalization matrices not loaded, skipped" << ENDL;
        DoNorm = false;
    }
    if(GatePx.first ==0.) {
        ERR_MESS << "You need first to define the gates using: SetGate(Px,w_Px,Py,w_Py)" << ENDL;
        return false;
    }

    if(!DoStdBgd && DoAdvBgd) {
        WARN_MESS << "Cannot process advanced background without standard, bgd skipped" << ENDL;
        DoAdvBgd = false;
    }

    // Check and adapt if necessary the widths of the gates to correspond to the closests bins

    vector<TBox *> boxes{boxPxPy,boxPxBy,boxBxPy,boxBxBy,boxBxB2y,boxB2xBy,boxB2xB2y};
    TH1 *hist_test = vListOfMatrices.front().begin()->second;
    for(auto box: boxes) {
        box->SetX1(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetXaxis()->FindBin(box->GetX1())));
        box->SetX2(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetXaxis()->FindBin(box->GetX2())));
        box->SetY1(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetYaxis()->FindBin(box->GetY1())));
        box->SetY2(hist_test->GetXaxis()->GetBinLowEdge(hist_test->GetYaxis()->FindBin(box->GetY2())));
    }

    if(!DoAdvBgd) {
        boxBxB2y->SetLineWidth(0);boxBxB2y->SetFillStyle(0);
        boxB2xBy->SetLineWidth(0);boxB2xBy->SetFillStyle(0);
        boxB2xB2y->SetLineWidth(0);boxB2xB2y->SetFillStyle(0);
    }
    else {
        boxBxB2y->SetFillStyle(3001);boxBxB2y->SetLineWidth(3);
        boxB2xBy->SetFillStyle(3001);boxB2xBy->SetLineWidth(3);
        boxB2xB2y->SetFillStyle(3001);boxB2xB2y->SetLineWidth(3);
    }

    fCanvas->GetPad(1)->Modified();
    fCanvas->GetPad(1)->Update();

    int colors[3] = {kRed,kBlue,kGreen};

    for(int imode=0 ; imode<3 ; imode++) {

        if(imode==0 && (fCurrentMode[0] != '1')) continue;
        if(imode==1 && (fCurrentMode[1] != '1')) continue;
        if(imode==2 && (fCurrentMode[2] != '1')) continue;

        delete gAngCorr.at(imode);
        delete gAngCorrNorm.at(imode);
        delete gAngCorrCorr.at(imode);

        gAngCorr.at(imode) = nullptr;
        gAngCorrNorm.at(imode) = nullptr;
        gAngCorrCorr.at(imode) = nullptr;

        gAngCorr.at(imode) = new TGraphErrors;
        gAngCorr.at(imode)->SetNameTitle(Form("Angular_Correlation_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode),Form("Angular_Correlation_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode));
        gAngCorr.at(imode)->SetMarkerColor(colors[imode]);

        gAngCorr.at(imode)->SetMarkerStyle(20);
        gAngCorr.at(imode)->GetXaxis()->SetTitle(Form("%s: Angle (degree)",DetTypesString[imode].Data()));
        gAngCorr.at(imode)->GetYaxis()->SetTitle("Counts");
        gAngCorr.at(imode)->GetXaxis()->SetTitleSize(0.05);
        gAngCorr.at(imode)->GetXaxis()->SetLabelSize(0.05);
        gAngCorr.at(imode)->GetYaxis()->SetTitleSize(0.05);
        gAngCorr.at(imode)->GetYaxis()->SetLabelSize(0.05);
        gAngCorr.at(imode)->GetXaxis()->CenterTitle();
        gAngCorr.at(imode)->GetYaxis()->CenterTitle();

        if(DoNorm) {
            gAngCorrNorm.at(imode) = new TGraphErrors;
            gAngCorrNorm.at(imode)->SetNameTitle(Form("Angular_Correlation_Normalization_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode),Form("Angular_Correlation_Normalization_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode));
            gAngCorrNorm.at(imode)->SetMarkerColor(colors[imode]);
            gAngCorrNorm.at(imode)->SetMarkerStyle(20);
            gAngCorrNorm.at(imode)->GetXaxis()->SetTitle(Form("%s: Angle (degree)",DetTypesString[imode].Data()));
            gAngCorrNorm.at(imode)->GetYaxis()->SetTitle("Counts");
            gAngCorrNorm.at(imode)->GetXaxis()->SetTitleSize(0.05);
            gAngCorrNorm.at(imode)->GetXaxis()->SetLabelSize(0.05);
            gAngCorrNorm.at(imode)->GetYaxis()->SetTitleSize(0.05);
            gAngCorrNorm.at(imode)->GetYaxis()->SetLabelSize(0.05);
            gAngCorrNorm.at(imode)->GetXaxis()->CenterTitle();
            gAngCorrNorm.at(imode)->GetYaxis()->CenterTitle();

            gAngCorrCorr.at(imode) = new TGraphErrors;
            gAngCorrCorr.at(imode)->SetNameTitle(Form("Angular_Correlation_Normalized_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode),Form("Angular_Correlation_Normalized_%.1f_%.1f_%d",GatePx.first,GatePy.first,imode));
            gAngCorrCorr.at(imode)->SetMarkerColor(colors[imode]);
            gAngCorrCorr.at(imode)->SetMarkerStyle(20);
            gAngCorrCorr.at(imode)->GetXaxis()->SetTitle(Form("%s: Angle (degree)",DetTypesString[imode].Data()));
            gAngCorrCorr.at(imode)->GetYaxis()->SetTitle("Normalized Counts");
            gAngCorrCorr.at(imode)->GetXaxis()->SetTitleSize(0.05);
            gAngCorrCorr.at(imode)->GetXaxis()->SetLabelSize(0.05);
            gAngCorrCorr.at(imode)->GetYaxis()->SetTitleSize(0.05);
            gAngCorrCorr.at(imode)->GetYaxis()->SetLabelSize(0.05);
            gAngCorrCorr.at(imode)->GetYaxis()->SetTitleOffset(1.5);
            gAngCorrCorr.at(imode)->GetXaxis()->CenterTitle();
            gAngCorrCorr.at(imode)->GetYaxis()->CenterTitle();
        }

        Double_t WPx = (boxPxPy->GetX2()-boxPxPy->GetX1());
        Double_t WPy = (boxPxPy->GetY2()-boxPxPy->GetY1());

        if(Print) {
            INFO_MESS << "Projection without normalization, mode " << imode << ENDL;
        }
        for(auto i : vListOfMatrices.at(imode)) {

            Int_t Theta = i.first;
            TH2 *hist = i.second;

            Double_t PxPy = 0.;
            Double_t PxBy = 0.;
            Double_t BxPy = 0.;
            Double_t BxBy = 0.;
            Double_t BxB2y = 0.;
            Double_t B2xBy = 0.;
            Double_t B2xB2y = 0.;

            Double_t Total = 0.;
            Double_t Total_Err = 0.;

            PxPy = hist->Integral(hist->GetXaxis()->FindBin(boxPxPy->GetX1()),hist->GetXaxis()->FindBin(boxPxPy->GetX2())-1,hist->GetYaxis()->FindBin(boxPxPy->GetY1()),hist->GetYaxis()->FindBin(boxPxPy->GetY2())-1);

            if(DoStdBgd && boxPxBy) {
                PxBy = hist->Integral(hist->GetXaxis()->FindBin(boxPxBy->GetX1()),hist->GetXaxis()->FindBin(boxPxBy->GetX2())-1,hist->GetYaxis()->FindBin(boxPxBy->GetY1()),hist->GetYaxis()->FindBin(boxPxBy->GetY2())-1);
                PxBy /= ((boxPxBy->GetX2()-boxPxBy->GetX1())/WPx*(boxPxBy->GetY2()-boxPxBy->GetY1())/WPy);
                BxPy = hist->Integral(hist->GetXaxis()->FindBin(boxBxPy->GetX1()),hist->GetXaxis()->FindBin(boxBxPy->GetX2())-1,hist->GetYaxis()->FindBin(boxBxPy->GetY1()),hist->GetYaxis()->FindBin(boxBxPy->GetY2())-1);
                BxPy /= ((boxBxPy->GetX2()-boxBxPy->GetX1())/WPx*(boxBxPy->GetY2()-boxBxPy->GetY1())/WPy);
                BxBy = hist->Integral(hist->GetXaxis()->FindBin(boxBxBy->GetX1()),hist->GetXaxis()->FindBin(boxBxBy->GetX2())-1,hist->GetYaxis()->FindBin(boxBxBy->GetY1()),hist->GetYaxis()->FindBin(boxBxBy->GetY2())-1);
                BxBy /= ((boxBxBy->GetX2()-boxBxBy->GetX1())/WPx*(boxBxBy->GetY2()-boxBxBy->GetY1())/WPy);
            }
            if(DoStdBgd && DoAdvBgd && boxBxB2y) {
                BxB2y = hist->Integral(hist->GetXaxis()->FindBin(boxBxB2y->GetX1()),hist->GetXaxis()->FindBin(boxBxB2y->GetX2())-1,hist->GetYaxis()->FindBin(boxBxB2y->GetY1()),hist->GetYaxis()->FindBin(boxBxB2y->GetY2())-1);
                BxB2y /= ((boxBxB2y->GetX2()-boxBxB2y->GetX1())/WPx*(boxBxB2y->GetY2()-boxBxB2y->GetY1())/WPy);
                B2xBy = hist->Integral(hist->GetXaxis()->FindBin(boxB2xBy->GetX1()),hist->GetXaxis()->FindBin(boxB2xBy->GetX2())-1,hist->GetYaxis()->FindBin(boxB2xBy->GetY1()),hist->GetYaxis()->FindBin(boxB2xBy->GetY2())-1);
                B2xBy /= ((boxB2xBy->GetX2()-boxB2xBy->GetX1())/WPx*(boxB2xBy->GetY2()-boxB2xBy->GetY1())/WPy);
                B2xB2y = hist->Integral(hist->GetXaxis()->FindBin(boxB2xB2y->GetX1()),hist->GetXaxis()->FindBin(boxB2xB2y->GetX2())-1,hist->GetYaxis()->FindBin(boxB2xB2y->GetY1()),hist->GetYaxis()->FindBin(boxB2xB2y->GetY2())-1);
                B2xB2y /= ((boxB2xB2y->GetX2()-boxB2xB2y->GetX1())/WPx*(boxB2xB2y->GetY2()-boxB2xB2y->GetY1())/WPy);
            }

            Total = PxPy - PxBy - BxPy + BxBy + 0.5*( -BxB2y -B2xBy + 2.*B2xB2y);
            Total_Err = 2*sqrt(PxPy + PxBy + BxPy + BxBy + 0.5*0.5*(BxB2y + B2xBy + B2xB2y));

            gAngCorr.at(imode)->SetPoint(gAngCorr.at(imode)->GetN(),Theta,Total);
            gAngCorr.at(imode)->SetPointError(gAngCorr.at(imode)->GetN()-1,2,Total_Err);
        }

        //        // Scaling
        //        Double_t Mean = 0.;
        //        for(int i=0; i<gAngCorr.at(imode)->GetN() ; i++) Mean += gAngCorr.at(imode)->GetY()[i];
        //        Mean /= (float)gAngCorr.at(imode)->GetN();
        //        for(int i=0; i<gAngCorr.at(imode)->GetN() ; i++) {
        //            gAngCorr.at(imode)->SetPoint(i,gAngCorr.at(imode)->GetX()[i],gAngCorr.at(imode)->GetY()[i]/Mean);
        //            gAngCorr.at(imode)->SetPointError(i,gAngCorr.at(imode)->GetEX()[i],gAngCorr.at(imode)->GetEY()[i]/Mean);
        //        }

        if(NormTo90) {
            // Scaling to 90°
            Double_t Mean = gAngCorr.at(imode)->Eval(90);
            for(int i=0; i<gAngCorr.at(imode)->GetN() ; i++) {
                gAngCorr.at(imode)->SetPoint(i,gAngCorr.at(imode)->GetX()[i],gAngCorr.at(imode)->GetY()[i]/Mean);
                gAngCorr.at(imode)->SetPointError(i,gAngCorr.at(imode)->GetEX()[i],gAngCorr.at(imode)->GetEY()[i]/Mean);
            }
        }

        if(Print) {
            for(int i=0; i<gAngCorr.at(imode)->GetN() ; i++) {
                INFO_MESS << " -> Theta = " << gAngCorr.at(imode)->GetX()[i] << ", N = " << gAngCorr.at(imode)->GetY()[i] << ", Err = " << gAngCorr.at(imode)->GetEY()[i] << ENDL;
            }
        }

        gAngCorr.at(imode)->GetXaxis()->SetRangeUser(0,180);

        if(DoNorm == 1) {
            if(Print) {
                INFO_MESS << "Projection of the normalization matrices, mode " << imode << ENDL;
            }
            for(auto i : vListOfNormMatrices.at(imode)) {
                Int_t Theta = i.first;
                TH2 *hist = i.second;

                Double_t PxPy = 0.;
                Double_t PxBy = 0.;
                Double_t BxPy = 0.;
                Double_t BxBy = 0.;
                Double_t BxB2y = 0.;
                Double_t B2xBy = 0.;
                Double_t B2xB2y = 0.;

                Double_t Total = 0.;
                Double_t Total_Err = 0.;

                PxPy = hist->Integral(hist->GetXaxis()->FindBin(boxPxPy->GetX1()),hist->GetXaxis()->FindBin(boxPxPy->GetX2())-1,hist->GetYaxis()->FindBin(boxPxPy->GetY1()),hist->GetYaxis()->FindBin(boxPxPy->GetY2())-1);

                if(DoStdBgd && boxPxBy) {
                    PxBy = hist->Integral(hist->GetXaxis()->FindBin(boxPxBy->GetX1()),hist->GetXaxis()->FindBin(boxPxBy->GetX2())-1,hist->GetYaxis()->FindBin(boxPxBy->GetY1()),hist->GetYaxis()->FindBin(boxPxBy->GetY2())-1);
                    PxBy /= ((boxPxBy->GetX2()-boxPxBy->GetX1())/WPx*(boxPxBy->GetY2()-boxPxBy->GetY1())/WPy);

                    BxPy = hist->Integral(hist->GetXaxis()->FindBin(boxBxPy->GetX1()),hist->GetXaxis()->FindBin(boxBxPy->GetX2())-1,hist->GetYaxis()->FindBin(boxBxPy->GetY1()),hist->GetYaxis()->FindBin(boxBxPy->GetY2())-1);
                    BxPy /= ((boxBxPy->GetX2()-boxBxPy->GetX1())/WPx*(boxBxPy->GetY2()-boxBxPy->GetY1())/WPy);

                    BxBy = hist->Integral(hist->GetXaxis()->FindBin(boxBxBy->GetX1()),hist->GetXaxis()->FindBin(boxBxBy->GetX2())-1,hist->GetYaxis()->FindBin(boxBxBy->GetY1()),hist->GetYaxis()->FindBin(boxBxBy->GetY2())-1);
                    BxBy /= ((boxBxBy->GetX2()-boxBxBy->GetX1())/WPx*(boxBxBy->GetY2()-boxBxBy->GetY1())/WPy);
                }
                if(DoStdBgd && DoAdvBgd && boxBxB2y) {
                    BxB2y = hist->Integral(hist->GetXaxis()->FindBin(boxBxB2y->GetX1()),hist->GetXaxis()->FindBin(boxBxB2y->GetX2())-1,hist->GetYaxis()->FindBin(boxBxB2y->GetY1()),hist->GetYaxis()->FindBin(boxBxB2y->GetY2())-1);
                    BxB2y /= ((boxBxB2y->GetX2()-boxBxB2y->GetX1())/WPx*(boxBxB2y->GetY2()-boxBxB2y->GetY1())/WPy);

                    B2xBy = hist->Integral(hist->GetXaxis()->FindBin(boxB2xBy->GetX1()),hist->GetXaxis()->FindBin(boxB2xBy->GetX2())-1,hist->GetYaxis()->FindBin(boxB2xBy->GetY1()),hist->GetYaxis()->FindBin(boxB2xBy->GetY2())-1);
                    B2xBy /= ((boxB2xBy->GetX2()-boxB2xBy->GetX1())/WPx*(boxB2xBy->GetY2()-boxB2xBy->GetY1())/WPy);

                    B2xB2y = hist->Integral(hist->GetXaxis()->FindBin(boxB2xB2y->GetX1()),hist->GetXaxis()->FindBin(boxB2xB2y->GetX2())-1,hist->GetYaxis()->FindBin(boxB2xB2y->GetY1()),hist->GetYaxis()->FindBin(boxB2xB2y->GetY2())-1);
                    B2xB2y /= ((boxB2xB2y->GetX2()-boxB2xB2y->GetX1())/WPx*(boxB2xB2y->GetY2()-boxB2xB2y->GetY1())/WPy);
                }

                Total = PxPy /*- PxBy - BxPy + BxBy + 0.5*( -BxB2y -B2xBy + 2.*B2xB2y)*/;
                Total_Err = 2*sqrt(PxPy /*+ PxBy + BxPy + BxBy + 0.5*0.5*(BxB2y + B2xBy + B2xB2y)*/);

                gAngCorrNorm.at(imode)->SetPoint(gAngCorrNorm.at(imode)->GetN(),Theta,Total);
                gAngCorrNorm.at(imode)->SetPointError(gAngCorrNorm.at(imode)->GetN()-1,2,Total_Err);
            }

            //            // Scaling to mean value
            //            Double_t Mean = 0.;
            //            for(int i=0; i<gAngCorrNorm.at(imode)->GetN() ; i++) Mean += gAngCorrNorm.at(imode)->GetY()[i];
            //            Mean /= (float)gAngCorrNorm.at(imode)->GetN();
            //            for(int i=0; i<gAngCorrNorm.at(imode)->GetN() ; i++) {
            //                gAngCorrNorm.at(imode)->SetPoint(i,gAngCorrNorm.at(imode)->GetX()[i],gAngCorrNorm.at(imode)->GetY()[i]/Mean);
            //                gAngCorrNorm.at(imode)->SetPointError(i,gAngCorrNorm.at(imode)->GetEX()[i],gAngCorrNorm.at(imode)->GetEY()[i]/Mean);
            //            }

            if(NormTo90) {
                // Scaling to 90°
                Double_t Mean = gAngCorrNorm.at(imode)->Eval(90);
                for(int i=0; i<gAngCorrNorm.at(imode)->GetN() ; i++) {
                    gAngCorrNorm.at(imode)->SetPoint(i,gAngCorrNorm.at(imode)->GetX()[i],gAngCorrNorm.at(imode)->GetY()[i]/Mean);
                    gAngCorrNorm.at(imode)->SetPointError(i,gAngCorrNorm.at(imode)->GetEX()[i],gAngCorrNorm.at(imode)->GetEY()[i]/Mean);
                }
            }

            if(Print) {
                for(int i=0; i<gAngCorrNorm.at(imode)->GetN() ; i++) {
                    INFO_MESS << " -> Theta = " << gAngCorrNorm.at(imode)->GetX()[i] << ", N = " << gAngCorrNorm.at(imode)->GetY()[i] << ", Err = " << gAngCorrNorm.at(imode)->GetEY()[i] << ENDL;
                }
            }

            for(int i=0; i<gAngCorr.at(imode)->GetN() ; i++) {
                gAngCorrCorr.at(imode)->SetPoint(i,gAngCorr.at(imode)->GetX()[i],gAngCorr.at(imode)->GetY()[i]/gAngCorrNorm.at(imode)->GetY()[i]);
                gAngCorrCorr.at(imode)->SetPointError(i,gAngCorr.at(imode)->GetEX()[i], gAngCorrCorr.at(imode)->GetY()[i]*sqrt((gAngCorr.at(imode)->GetEY()[i]/gAngCorr.at(imode)->GetY()[i])*(gAngCorr.at(imode)->GetEY()[i]/gAngCorr.at(imode)->GetY()[i]) +(gAngCorrNorm.at(imode)->GetEY()[i]/gAngCorrNorm.at(imode)->GetY()[i])*(gAngCorrNorm.at(imode)->GetEY()[i]/gAngCorrNorm.at(imode)->GetY()[i])));
            }

            gAngCorrCorr.at(imode)->GetXaxis()->SetRangeUser(0,180);
        }
    }

    // adapt all axis
    double ymin=1e6;
    double ymax=-1e6;
    for(int imode=0 ; imode<3 ; imode++) {
        if(imode==0 && (fCurrentMode[0] != '1')) continue;
        if(imode==1 && (fCurrentMode[1] != '1')) continue;
        if(imode==2 && (fCurrentMode[2] != '1')) continue;

        if(gAngCorrCorr.at(imode)->GetYaxis()->GetXmax()>ymax) ymax = gAngCorrCorr.at(imode)->GetYaxis()->GetXmax();
        if(gAngCorrCorr.at(imode)->GetYaxis()->GetXmin()<ymin) ymin = gAngCorrCorr.at(imode)->GetYaxis()->GetXmin();
    }
    for(int imode=0 ; imode<3 ; imode++) {
        if(imode==0 && (fCurrentMode[0] != '1')) continue;
        if(imode==1 && (fCurrentMode[1] != '1')) continue;
        if(imode==2 && (fCurrentMode[2] != '1')) continue;

        gAngCorrCorr.at(imode)->GetYaxis()->SetRangeUser(ymin,ymax);
    }

    return true;
}

double AngCorrMinimizer(const double *xx)
{
    Double_t Chi2 = 0;

    Double_t A0  = xx[0];
    Double_t A2  = xx[1];
    Double_t A4  = xx[2];

    Double_t Q2[3];
    Double_t Q4[3];
    for(int i=0 ; i<3 ; i++) {
        Q2[i] = xx[3+2*i];
        Q4[i] = xx[3+2*i+1];
    }
    for(int imode=0 ; imode<3 ; imode++) {
        if(imode==0 && (fCurrentMode[0] != '1')) continue;
        if(imode==1 && (fCurrentMode[1] != '1')) continue;
        if(imode==2 && (fCurrentMode[2] != '1')) continue;

        for(int i=0 ; i<gAngCorrCorr.at(imode)->GetN() ; i++) {
            Double_t Theta = gAngCorrCorr.at(imode)->GetX()[i];
            //        Double_t eTheta = gAngCorrCorr->GetEX()[i];

            Double_t W = gAngCorrCorr.at(imode)->GetY()[i];
            Double_t eW = gAngCorrCorr.at(imode)->GetEY()[i];

            Double_t f_Theta = A0*(1+A2*Q2[imode]*ROOT::Math::legendre(2,cos(Theta*TMath::DegToRad())) + A4*Q4[imode]*ROOT::Math::legendre(4,cos(Theta*TMath::DegToRad())));

            Chi2 += (W-f_Theta)*(W-f_Theta)/(eW*eW);
        }
    }

    return Chi2;
}

void Minimize(bool fixQ2 = true, bool fixQ4 = true, bool fixA=false, double a2=0.5, double a4=0.1)
{
    // create minimizer giving a name and a name (optionally) for the specific
    // algorithm
    // possible choices are:
    //     minName                  algoName
    // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
    //  Minuit2                     Fumili2
    //  Fumili
    //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
    //                              BFGS2, SteepestDescent
    //  GSLMultiFit
    //   GSLSimAn
    //   Genetic
    delete fMinimum;
    fMinimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerance , etc...
    fMinimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    fMinimum->SetMaxIterations(10000);  // for GSL
    fMinimum->SetTolerance(0.001);
    fMinimum->SetPrintLevel(1);

    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(&AngCorrMinimizer,9);

    double step[9]/* = {0.01,0.001,0.001,0.01,0.}*/;
    // define steps
    // A0 A1 A2
    step[0] = 0.01;
    step[1] = 0.001;
    step[2] = 0.001;
    // Qvalues
    for(int i=0 ; i<3 ; i++) {
        step[3+2*i]   = 0.01;
        step[3+2*i+1] = 0.01;
    }

    // starting point
    double variable[9];
    // A0 A1 A2
    variable[0] = 1.;
    variable[1] = a2;
    variable[2] = a4;
    // Qvalues
    for(int i=0 ; i<3 ; i++) {
        variable[3+2*i]   = fQValues.at(i).first.first;
        variable[3+2*i+1] = fQValues.at(i).second.first;
    }

    fMinimum->SetFunction(f);

    // Set the free variables to be minimized !
    fMinimum->SetVariable(0,"A0",variable[0], step[0]);
    fMinimum->SetVariable(1,"A2",variable[1], step[1]);
    fMinimum->SetVariable(2,"A4",variable[2], step[2]);
    for(int i=0 ; i<3 ; i++) {
        fMinimum->SetLimitedVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),variable[3+2*i],   step[3+2*i],0.,1.);
        fMinimum->SetLimitedVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),variable[3+2*i+1], step[3+2*i+1],0.,1.);
    }

    if(fixQ2 || fixQ4) {
        for(int i=0 ; i<3 ; i++) {
            if(fixQ2) fMinimum->FixVariable(3+2*i);
            if(fixQ4) fMinimum->FixVariable(3+2*i+1);
        }
    }
    if(fixA) {
        fMinimum->FixVariable(1);
        fMinimum->FixVariable(2);
    }

    //check the mode of the fit
    for(int i=0 ; i<3 ; i++) {
        if(fCurrentMode[i] == '0') {
            fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),0.,0.);
            fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),0.,0.);
            fMinimum->FixVariable(3+2*i);
            fMinimum->FixVariable(3+2*i+1);
        }
    }

    // do the minimization
    fMinimum->SetPrintLevel(0);
    fMinimum->Minimize();


    if(fixQ2 && fixQ4) {

        // extract the fluctuations on A2 and A4 using the errors on Q2i and Q4i
        double minA0 =  1e12;
        double maxA0 =  -1e12;
        double minA2 =  1e12;
        double maxA2 =  -1e12;
        double minA4 =  1e12;
        double maxA4 =  -1e12;

        if(fUseMonteCarlo) {
            TStatistic stat_A0;
            TStatistic stat_A2;
            TStatistic stat_A4;

            TString MCString = Form("%d",fNMonteCarlo);
            for(int itest=0 ; itest<fNMonteCarlo ; itest++) {
                cout<< left << "Monte-Carlo estimation of errors: "<< setw(MCString.Length()) << itest<<"/"<<fNMonteCarlo<< " :" << ((int)((double)itest)/(double(fNMonteCarlo))*100.) << "\%:" << "\r";
                for(int i=0 ; i<3 ; i++) {
                    if(fCurrentMode[i] == '1') {
                        double _Q2 = gRandom->Uniform(fQValues.at(i).first.first-fQValues.at(i).first.second,fQValues.at(i).first.first+fQValues.at(i).first.second);
                        double _Q4 = gRandom->Uniform(fQValues.at(i).second.first-fQValues.at(i).second.second,fQValues.at(i).second.first+fQValues.at(i).second.second);
                        if(_Q2>1.) _Q2=1.;
                        if(_Q4>1.) _Q4=1.;
                        fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),_Q2, step[3+2*i]);
                        fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),_Q4, step[3+2*i+1]);
                        fMinimum->FixVariable(3+2*i);
                        fMinimum->FixVariable(3+2*i+1);
                    }
                }
                fMinimum->Minimize();
                stat_A0.Fill(fMinimum->X()[0]);
                stat_A2.Fill(fMinimum->X()[1]);
                stat_A4.Fill(fMinimum->X()[2]);

                if((fMinimum->X()[0]-fMinimum->Errors()[0])<minA0) minA0 = (fMinimum->X()[0]-fMinimum->Errors()[0]);
                if((fMinimum->X()[0]+fMinimum->Errors()[0])>maxA0) maxA0 = (fMinimum->X()[0]+fMinimum->Errors()[0]);
                if((fMinimum->X()[1]-fMinimum->Errors()[1])<minA2) minA2 = (fMinimum->X()[1]-fMinimum->Errors()[1]);
                if((fMinimum->X()[1]+fMinimum->Errors()[1])>maxA2) maxA2 = (fMinimum->X()[1]+fMinimum->Errors()[1]);
                if((fMinimum->X()[2]-fMinimum->Errors()[2])<minA4) minA4 = (fMinimum->X()[2]-fMinimum->Errors()[2]);
                if((fMinimum->X()[2]+fMinimum->Errors()[2])>maxA4) maxA4 = (fMinimum->X()[2]+fMinimum->Errors()[2]);
            }
            //            maxA0 = stat_A0.GetMax();
            //            minA0 = stat_A0.GetMin();
            //            maxA2 = stat_A2.GetMax();
            //            minA2 = stat_A2.GetMin();
            //            maxA4 = stat_A4.GetMax();
            //            minA4 = stat_A4.GetMin();
        }
        else {
            // Q2 analysis
            for(int i=0 ; i<3 ; i++) {
                if(fCurrentMode[i] == '1') {
                    double _Q2 = fQValues.at(i).first.first  - fQValues.at(i).first.second;
                    double _Q4 = fQValues.at(i).second.first - fQValues.at(i).second.second;
                    if(_Q2>1.) _Q2=1.;
                    if(_Q4>1.) _Q4=1.;
                    fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),_Q2, step[3+2*i]);
                    fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),_Q4, step[3+2*i+1]);
                    fMinimum->FixVariable(3+2*i);
                    fMinimum->FixVariable(3+2*i+1);
                }
            }
            fMinimum->Minimize();
            if((fMinimum->X()[0]-fMinimum->Errors()[0])<minA0) minA0 = (fMinimum->X()[0]-fMinimum->Errors()[0]);
            if((fMinimum->X()[0]+fMinimum->Errors()[0])>maxA0) maxA0 = (fMinimum->X()[0]+fMinimum->Errors()[0]);
            if((fMinimum->X()[1]-fMinimum->Errors()[1])<minA2) minA2 = (fMinimum->X()[1]-fMinimum->Errors()[1]);
            if((fMinimum->X()[1]+fMinimum->Errors()[1])>maxA2) maxA2 = (fMinimum->X()[1]+fMinimum->Errors()[1]);
            if((fMinimum->X()[2]-fMinimum->Errors()[2])<minA4) minA4 = (fMinimum->X()[2]-fMinimum->Errors()[2]);
            if((fMinimum->X()[2]+fMinimum->Errors()[2])>maxA4) maxA4 = (fMinimum->X()[2]+fMinimum->Errors()[2]);

            for(int i=0 ; i<3 ; i++) {
                if(fCurrentMode[i] == '1') {
                    double _Q2 = fQValues.at(i).first.first  + fQValues.at(i).first.second;
                    double _Q4 = fQValues.at(i).second.first + fQValues.at(i).second.second;
                    if(_Q2>1.) _Q2=1.;
                    if(_Q4>1.) _Q4=1.;
                    fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),_Q2, step[3+2*i]);
                    fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),_Q4, step[3+2*i+1]);
                    fMinimum->FixVariable(3+2*i);
                    fMinimum->FixVariable(3+2*i+1);
                }
            }
            fMinimum->Minimize();
            if((fMinimum->X()[0]-fMinimum->Errors()[0])<minA0) minA0 = (fMinimum->X()[0]-fMinimum->Errors()[0]);
            if((fMinimum->X()[0]+fMinimum->Errors()[0])>maxA0) maxA0 = (fMinimum->X()[0]+fMinimum->Errors()[0]);
            if((fMinimum->X()[1]-fMinimum->Errors()[1])<minA2) minA2 = (fMinimum->X()[1]-fMinimum->Errors()[1]);
            if((fMinimum->X()[1]+fMinimum->Errors()[1])>maxA2) maxA2 = (fMinimum->X()[1]+fMinimum->Errors()[1]);
            if((fMinimum->X()[2]-fMinimum->Errors()[2])<minA4) minA4 = (fMinimum->X()[2]-fMinimum->Errors()[2]);
            if((fMinimum->X()[2]+fMinimum->Errors()[2])>maxA4) maxA4 = (fMinimum->X()[2]+fMinimum->Errors()[2]);

            for(int i=0 ; i<3 ; i++) {
                if(fCurrentMode[i] == '1') {
                    double _Q2 = fQValues.at(i).first.first  + fQValues.at(i).first.second;
                    double _Q4 = fQValues.at(i).second.first - fQValues.at(i).second.second;
                    if(_Q2>1.) _Q2=1.;
                    if(_Q4>1.) _Q4=1.;
                    fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),_Q2, step[3+2*i]);
                    fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),_Q4, step[3+2*i+1]);
                    fMinimum->FixVariable(3+2*i);
                    fMinimum->FixVariable(3+2*i+1);
                }
            }
            fMinimum->Minimize();
            if((fMinimum->X()[0]-fMinimum->Errors()[0])<minA0) minA0 = (fMinimum->X()[0]-fMinimum->Errors()[0]);
            if((fMinimum->X()[0]+fMinimum->Errors()[0])>maxA0) maxA0 = (fMinimum->X()[0]+fMinimum->Errors()[0]);
            if((fMinimum->X()[1]-fMinimum->Errors()[1])<minA2) minA2 = (fMinimum->X()[1]-fMinimum->Errors()[1]);
            if((fMinimum->X()[1]+fMinimum->Errors()[1])>maxA2) maxA2 = (fMinimum->X()[1]+fMinimum->Errors()[1]);
            if((fMinimum->X()[2]-fMinimum->Errors()[2])<minA4) minA4 = (fMinimum->X()[2]-fMinimum->Errors()[2]);
            if((fMinimum->X()[2]+fMinimum->Errors()[2])>maxA4) maxA4 = (fMinimum->X()[2]+fMinimum->Errors()[2]);

            for(int i=0 ; i<3 ; i++) {
                if(fCurrentMode[i] == '1') {
                    double _Q2 = fQValues.at(i).first.first  - fQValues.at(i).first.second;
                    double _Q4 = fQValues.at(i).second.first + fQValues.at(i).second.second;
                    if(_Q2>1.) _Q2=1.;
                    if(_Q4>1.) _Q4=1.;
                    fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),_Q2, step[3+2*i]);
                    fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),_Q4, step[3+2*i+1]);
                    fMinimum->FixVariable(3+2*i);
                    fMinimum->FixVariable(3+2*i+1);
                }
            }
            fMinimum->Minimize();
            if((fMinimum->X()[0]-fMinimum->Errors()[0])<minA0) minA0 = (fMinimum->X()[0]-fMinimum->Errors()[0]);
            if((fMinimum->X()[0]+fMinimum->Errors()[0])>maxA0) maxA0 = (fMinimum->X()[0]+fMinimum->Errors()[0]);
            if((fMinimum->X()[1]-fMinimum->Errors()[1])<minA2) minA2 = (fMinimum->X()[1]-fMinimum->Errors()[1]);
            if((fMinimum->X()[1]+fMinimum->Errors()[1])>maxA2) maxA2 = (fMinimum->X()[1]+fMinimum->Errors()[1]);
            if((fMinimum->X()[2]-fMinimum->Errors()[2])<minA4) minA4 = (fMinimum->X()[2]-fMinimum->Errors()[2]);
            if((fMinimum->X()[2]+fMinimum->Errors()[2])>maxA4) maxA4 = (fMinimum->X()[2]+fMinimum->Errors()[2]);
        }

        // set the real parameters
        for(int i=0 ; i<3 ; i++) {
            if(fCurrentMode[i] == '1') {
                // set Q2-Q2Err
                fMinimum->SetVariable(3+2*i,  Form("%s_Q2",DetTypesString[i].Data()),fQValues.at(i).first.first,   step[3+2*i]);
                fMinimum->SetVariable(3+2*i+1,Form("%s_Q4",DetTypesString[i].Data()),fQValues.at(i).second.first, step[3+2*i+1]);
                fMinimum->FixVariable(3+2*i);
                fMinimum->FixVariable(3+2*i+1);
            }
        }

        fMinimum->Minimize();

        fA0Exp = fMinimum->X()[0];
        fA0Exp_errlow  = fA0Exp - minA0;
        fA0Exp_errhigh = maxA0 - fA0Exp;

        fA2Exp = fMinimum->X()[1];
        fA2Exp_errlow  = fA2Exp - minA2;
        fA2Exp_errhigh = maxA2 - fA2Exp;

        fA4Exp = fMinimum->X()[2];
        fA4Exp_errlow  = fA4Exp - minA4;
        fA4Exp_errhigh = maxA4 - fA4Exp;
    }

    cout << "W(theta) = A0 x ( 1 + A2*Qi2*Pk2 + A4*Qi4*Pk4 ), using:" << endl;
    if(fixA) {
        cout << "=> A2 = " << a2 << endl;
        cout << "=> A4 = " << a4 << endl;
        for(int i=0 ; i<3 ; i++) {
            if(fCurrentMode[i] =='0') continue;
            if(fixQ2) cout << Form("=> %s Q2 = %.4f +- %.4f",DetTypesString[i].Data(),fQValues.at(i).first.first,fQValues.at(i).first.second) << endl;
            if(fixQ4) cout << Form("=> %s Q4 = %.4f +- %.4f",DetTypesString[i].Data(),fQValues.at(i).second.first,fQValues.at(i).second.second) << endl;
        }
        cout <<endl;
        cout << "Result: Chi2/ndf = " << fMinimum->MinValue()/((double)f.NDim())<<endl;
        cout << "             ndf = " << f.NDim() <<endl;
        cout << "=> A0 = "<< left << setw(10) << fMinimum->X()[0] << " +- " << fMinimum->Errors()[0] << endl;
        for(int i=0 ; i<3 ; i++) {
            if(fCurrentMode[i] =='0') continue;
            if(!fixQ2) cout << Form("=> %s Q2 = %.4f +- %.4f",DetTypesString[i].Data(),fMinimum->X()[3+2*i],fMinimum->Errors()[3+2*i]) << endl;
            if(!fixQ4) cout << Form("=> %s Q4 = %.4f +- %.4f",DetTypesString[i].Data(),fMinimum->X()[3+2*i+1],fMinimum->Errors()[3+2*i+1]) << endl;
        }
    }
    else {
        for(int i=0 ; i<3 ; i++) {
            if(fCurrentMode[i] =='0') continue;
            cout << Form("=> %s Q2 = %.4f +- %.4f",DetTypesString[i].Data(),fQValues.at(i).first.first,fQValues.at(i).first.second) << endl;
            cout << Form("=> %s Q4 = %.4f +- %.4f",DetTypesString[i].Data(),fQValues.at(i).second.first,fQValues.at(i).second.second) << endl;
        }
        cout <<endl;
        cout << "Result: Chi2/ndf = " << fMinimum->MinValue()/((double)f.NDim())<<endl;
        cout << "             ndf = " << f.NDim() <<endl;
        //        cout << "=> A0 = "<< left << setw(10) << fMinimum->X()[0] << " +- " << fMinimum->Errors()[0] << endl;
        //        cout << "=> A2 = "<< left << setw(10) << fMinimum->X()[1] << " +- " << fMinimum->Errors()[1] << endl;
        //        cout << "=> A4 = "<< left << setw(10) << fMinimum->X()[2] << " +- " << fMinimum->Errors()[2] << endl;

        cout << "=> A0 = "<< left << setw(10) << fMinimum->X()[0] << " [ - " << setw(10) << fA0Exp_errlow << " ; + " << setw(10) << fA0Exp_errhigh << " ]" << endl;
        cout << "=> A2 = "<< left << setw(10) << fMinimum->X()[1] << " [ - " << setw(10) << fA2Exp_errlow << " ; + " << setw(10) << fA2Exp_errhigh << " ]" << endl;
        cout << "=> A4 = "<< left << setw(10) << fMinimum->X()[2] << " [ - " << setw(10) << fA4Exp_errlow << " ; + " << setw(10) << fA4Exp_errhigh << " ]" << endl;
    }
}

void SetQValues()
{
    for(int i=0 ; i<3 ; i++) {
        if(fCurrentMode[i] =='0') continue;
        fQValues.at(i).first.first = fMinimum->X()[3+2*i];
        if(fMinimum->Errors()[3+2*i]>0.) fQValues.at(i).first.second = fMinimum->Errors()[3+2*i];
        fQValues.at(i).second.first = fMinimum->X()[3+2*i+1];
        if(fMinimum->Errors()[3+2*i+1]>0) fQValues.at(i).second.second = fMinimum->Errors()[3+2*i+1];
    }
    cout << left << setw(12) << "Det types  : " << setw(10) << "Q2" << setw(10) << "Q2 unc" << setw(10) << "Q4" << setw(10) << "Q4 unc" << endl;
    cout << left << setw(12) << "FIPPS-FIPPS: " << setw(10) << fQValues.at(kMode_FF).first.first << setw(10) << fQValues.at(kMode_FF).first.second << setw(10) << fQValues.at(kMode_FF).second.first << setw(10) << fQValues.at(kMode_FF).second.second << endl;
    cout << left << setw(12) << "FIPPS-IFIN : " << setw(10) << fQValues.at(kMode_FI).first.first << setw(10) << fQValues.at(kMode_FI).first.second << setw(10) << fQValues.at(kMode_FI).second.first << setw(10) << fQValues.at(kMode_FI).second.second << endl;
    cout << left << setw(12) << "IFIN-IFIN  : " << setw(10) << fQValues.at(kMode_II).first.first << setw(10) << fQValues.at(kMode_II).first.second << setw(10) << fQValues.at(kMode_II).second.first << setw(10) << fQValues.at(kMode_II).second.second << endl;
}

void SaveQValues()
{
    ofstream qval_file_o("qvalues.txt");
    qval_file_o << left << setw(12) << "Det types  : " << setw(10) << "Q2" << setw(10) << "Q2 unc" << setw(10) << "Q4" << setw(10) << "AQ unc" << endl;
    qval_file_o << left << setw(12) << "FIPPS-FIPPS: " << setw(10) << fQValues.at(kMode_FF).first.first << setw(10) << fQValues.at(kMode_FF).first.second << setw(10) << fQValues.at(kMode_FF).second.first << setw(10) << fQValues.at(kMode_FF).second.second << endl;
    qval_file_o << left << setw(12) << "FIPPS-IFIN : " << setw(10) << fQValues.at(kMode_FI).first.first << setw(10) << fQValues.at(kMode_FI).first.second << setw(10) << fQValues.at(kMode_FI).second.first << setw(10) << fQValues.at(kMode_FI).second.second << endl;
    qval_file_o << left << setw(12) << "IFIN-IFIN  : " << setw(10) << fQValues.at(kMode_II).first.first << setw(10) << fQValues.at(kMode_II).first.second << setw(10) << fQValues.at(kMode_II).second.first << setw(10) << fQValues.at(kMode_II).second.second << endl;
    qval_file_o.close();
}

void Project(bool DoStdBgd=1, bool DoAdvBgd=1, bool DoNorm=true, bool NormTo90=true, bool Print=false)
{
    bool ok = ProcessProject(DoNorm,DoStdBgd,DoAdvBgd,NormTo90,Print);
    if(!ok) return;
    if(DoNorm==false) {

        Int_t padpos=1;

        if(fCurrentMode[0] == '1') {
            fCanvas->cd(3)->cd(1)->cd(padpos);
            gAngCorr.at(kMode_FF)->Draw("ap");
            cout<<"FIPPS-FIPPS projection (no normalization):"<<endl;
            gAngCorr.at(kMode_FF)->Print();
            padpos++;
        }
        if(fCurrentMode[1] == '1') {
            fCanvas->cd(3)->cd(1)->cd(padpos);
            gAngCorr.at(kMode_FI)->Draw("ap");
            cout<<"FIPPS-IFIN projection (no normalization):"<<endl;
            gAngCorr.at(kMode_FI)->Print();
            padpos++;
        }
        if(fCurrentMode[2] == '1') {
            fCanvas->cd(3)->cd(1)->cd(padpos);
            gAngCorr.at(kMode_II)->Draw("ap");
            cout<<"IFIN-IFIN projection (no normalization):"<<endl;
            gAngCorr.at(kMode_II)->Print();
            padpos++;
        }
        return;
    }

    Minimize(true,true,false);

    for(int i=0 ; i<3 ; i++) {
        if(fCurrentMode[i] == '0') continue;

        delete fAngCorFuncFit.at(i);
        fAngCorFuncFit.at(i) = new TF1(Form("AngCorFuncFit_%d",i),AngCorrFuncion,0, 180,3);
        fAngCorFuncFit.at(i)->SetNpx(500);
        fAngCorFuncFit.at(i)->SetLineColor(gAngCorr.at(i)->GetMarkerColor());
        fAngCorFuncFit.at(i)->SetParameters(fMinimum->X()[0],fMinimum->X()[1]*fMinimum->X()[3+2*i],fMinimum->X()[2]*fMinimum->X()[3+2*i+1]);


        delete fAngCorFuncFit_Err.at(i);
        fAngCorFuncFit_Err.at(i) = new TGraphAsymmErrors;

        if(fUseMonteCarlo) {
            for(int ii=0 ; ii<180 ; ii++) {
                TStatistic stat;

                fAngCorFuncFit.at(i)->SetParameters(fMinimum->X()[0],fMinimum->X()[1]*fMinimum->X()[3+2*i],fMinimum->X()[2]*fMinimum->X()[3+2*i+1]);
                fAngCorFuncFit_Err.at(i)->SetPoint(ii,ii,fAngCorFuncFit.at(i)->Eval(ii));

                for(int test=0 ; test<fNMonteCarlo ; test++) {
                    double _A0 = gRandom->Uniform(fMinimum->X()[0]-fMinimum->Errors()[0],fMinimum->X()[0]+fMinimum->Errors()[0]);
                    double _A2 = gRandom->Uniform(fMinimum->X()[1]-fMinimum->Errors()[1],fMinimum->X()[1]+fMinimum->Errors()[1]);
                    double _A4 = gRandom->Uniform(fMinimum->X()[2]-fMinimum->Errors()[2],fMinimum->X()[2]+fMinimum->Errors()[2]);
                    double _Q2 = gRandom->Uniform(fQValues.at(i).first.first-fQValues.at(i).first.second,fQValues.at(i).first.first+fQValues.at(i).first.second);
                    double _Q4 = gRandom->Uniform(fQValues.at(i).second.first-fQValues.at(i).second.second,fQValues.at(i).second.first+fQValues.at(i).second.second);
                    if(_Q2>1.) _Q2=1.;
                    if(_Q4>1.) _Q4=1.;
                    fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                    stat.Fill(fAngCorFuncFit.at(i)->Eval(ii));
                }

                fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-stat.GetMin());
                fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,stat.GetMax()-fAngCorFuncFit_Err.at(i)->GetY()[ii]);
            }
        }
        else {
            for(int ii=0 ; ii<180 ; ii++) {
                fAngCorFuncFit.at(i)->SetParameters(fMinimum->X()[0],fMinimum->X()[1]*fMinimum->X()[3+2*i],fMinimum->X()[2]*fMinimum->X()[3+2*i+1]);
                fAngCorFuncFit_Err.at(i)->SetPoint(ii,ii,fAngCorFuncFit.at(i)->Eval(ii));

                double _A0 = fMinimum->X()[0]     + fMinimum->Errors()[0];
                double _A2 = fMinimum->X()[1]     - fMinimum->Errors()[1];
                double _Q2 = fMinimum->X()[3+2*i] + fQValues.at(i).first.second;
                double _A4 = fMinimum->X()[2]     - fMinimum->Errors()[0];
                double _Q4 = fMinimum->X()[3+2*i] + fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     + fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     + fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] - fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     + fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] - fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     + fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     - fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] + fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     + fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] - fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     + fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     + fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] - fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     - fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] + fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     - fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     - fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] + fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     - fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] + fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     - fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     + fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] - fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     + fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] - fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     - fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     - fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] + fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     + fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] - fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);

                _A0 = fMinimum->X()[0]     - fMinimum->Errors()[0];
                _A2 = fMinimum->X()[1]     + fMinimum->Errors()[1];
                _Q2 = fMinimum->X()[3+2*i] - fQValues.at(i).first.second;
                _A4 = fMinimum->X()[2]     - fMinimum->Errors()[0];
                _Q4 = fMinimum->X()[3+2*i] + fQValues.at(i).second.second;
                if(_Q2>1.) _Q2=1.;
                if(_Q4>1.) _Q4=1.;

                fAngCorFuncFit.at(i)->SetParameters(_A0,_A2*_Q2,_A4*_Q4);
                if(fAngCorFuncFit.at(i)->Eval(ii) < fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit_Err.at(i)->GetEYlow()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYlow(ii,fAngCorFuncFit_Err.at(i)->GetY()[ii]-fAngCorFuncFit.at(i)->Eval(ii));
                if(fAngCorFuncFit.at(i)->Eval(ii) > fAngCorFuncFit_Err.at(i)->GetY()[ii]+fAngCorFuncFit_Err.at(i)->GetEYhigh()[ii])
                    fAngCorFuncFit_Err.at(i)->SetPointEYhigh(ii,fAngCorFuncFit.at(i)->Eval(ii)-fAngCorFuncFit_Err.at(i)->GetY()[ii]);
            }
        }
        fAngCorFuncFit.at(i)->SetParameters(fMinimum->X()[0],fMinimum->X()[1]*fMinimum->X()[3+2*i],fMinimum->X()[2]*fMinimum->X()[3+2*i+1]);

        fAngCorFuncFit_Err.at(i)->SetFillStyle(3002);
        fAngCorFuncFit_Err.at(i)->SetFillColor(gAngCorr.at(i)->GetMarkerColor());
        fAngCorFuncFit_Err.at(i)->SetFillColorAlpha(gAngCorr.at(i)->GetMarkerColor(),0.5);
        fAngCorFuncFit_Err.at(i)->SetMarkerSize(0);
    }

    delete fAngCorrFunc;
    fAngCorrFunc = new TF1("AngCorrFunc",AngCorrFuncion,0, 180,3);
    fAngCorrFunc->SetNpx(500);
    fAngCorrFunc->SetLineColor(kBlack);
    fAngCorrFunc->SetParameters(fMinimum->X()[0],fMinimum->X()[1],fMinimum->X()[2]);
    fAngCorrFunc->SetParError(0,fMinimum->Errors()[0]);
    fAngCorrFunc->SetParError(1,fMinimum->Errors()[1]);
    fAngCorrFunc->SetParError(2,fMinimum->Errors()[2]);
    fAngCorrFunc->GetXaxis()->SetTitle("#theta (degree)");
    fAngCorrFunc->GetYaxis()->SetTitle("W(#theta)");
    fAngCorrFunc->GetYaxis()->SetTitleOffset(0.45);
    fAngCorrFunc->GetXaxis()->SetTitleSize(0.05);
    fAngCorrFunc->GetXaxis()->SetLabelSize(0.05);
    fAngCorrFunc->GetYaxis()->SetTitleSize(0.05);
    fAngCorrFunc->GetYaxis()->SetLabelSize(0.05);
    fAngCorrFunc->GetXaxis()->CenterTitle();
    fAngCorrFunc->GetYaxis()->CenterTitle();

    delete fAngCorrFuncErr;
    fAngCorrFuncErr = new TGraphAsymmErrors;

    if(fUseMonteCarlo) {
        for(int ii=0 ; ii<180 ; ii++) {
            TStatistic stat;

            fAngCorrFunc->SetParameters(fMinimum->X()[0],fMinimum->X()[1],fMinimum->X()[2]);
            fAngCorrFuncErr->SetPoint(ii,ii,fAngCorrFunc->Eval(ii));

            for(int test=0 ; test<fNMonteCarlo ; test++) {
                double _A0 = gRandom->Uniform(fA0Exp-fA0Exp_errlow,fA0Exp+fA0Exp_errhigh);
                double _A2 = gRandom->Uniform(fA2Exp-fA2Exp_errlow,fA2Exp+fA2Exp_errhigh);
                double _A4 = gRandom->Uniform(fA4Exp-fA4Exp_errlow,fA4Exp+fA4Exp_errhigh);

                fAngCorrFunc->SetParameters(_A0,_A2,_A4);
                stat.Fill(fAngCorrFunc->Eval(ii));
            }
            fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-stat.GetMin());
            fAngCorrFuncErr->SetPointEYhigh(ii,stat.GetMax()-fAngCorrFuncErr->GetY()[ii]);
        }
    }
    else {
        for(int ii=0 ; ii<180 ; ii++) {

            fAngCorrFunc->SetParameters(fMinimum->X()[0],fMinimum->X()[1],fMinimum->X()[2]);
            fAngCorrFuncErr->SetPoint(ii,ii,fAngCorrFunc->Eval(ii));

            double _A0 = fA0Exp + fA0Exp_errhigh;
            double _A2 = fA2Exp + fA2Exp_errhigh;
            double _A4 = fA4Exp + fA4Exp_errhigh;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);

            _A0 = fA0Exp - fA0Exp_errlow;
            _A2 = fA2Exp - fA2Exp_errlow;
            _A4 = fA4Exp - fA4Exp_errlow;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);

            _A0 = fA0Exp - fA0Exp_errlow;
            _A2 = fA2Exp + fA2Exp_errlow;
            _A4 = fA4Exp - fA4Exp_errlow;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);

            _A0 = fA0Exp - fA0Exp_errlow;
            _A2 = fA2Exp - fA2Exp_errlow;
            _A4 = fA4Exp + fA4Exp_errlow;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);

            _A0 = fA0Exp + fA0Exp_errlow;
            _A2 = fA2Exp + fA2Exp_errlow;
            _A4 = fA4Exp - fA4Exp_errlow;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);

            _A0 = fA0Exp + fA0Exp_errlow;
            _A2 = fA2Exp - fA2Exp_errlow;
            _A4 = fA4Exp + fA4Exp_errlow;

            fAngCorrFunc->SetParameters(_A0,_A2,_A4);

            if(fAngCorrFunc->Eval(ii) < fAngCorrFuncErr->GetY()[ii]-fAngCorrFuncErr->GetEYlow()[ii])
                fAngCorrFuncErr->SetPointEYlow(ii,fAngCorrFuncErr->GetY()[ii]-fAngCorrFunc->Eval(ii));
            if(fAngCorrFunc->Eval(ii) > fAngCorrFuncErr->GetY()[ii]+fAngCorrFuncErr->GetEYhigh()[ii])
                fAngCorrFuncErr->SetPointEYhigh(ii,fAngCorrFunc->Eval(ii)-fAngCorrFuncErr->GetY()[ii]);
        }
    }

    fAngCorrFuncErr->SetFillStyle(3002);
    fAngCorrFuncErr->SetFillColor(kBlack);
    fAngCorrFuncErr->SetFillColorAlpha(kBlack,0.5);
    fAngCorrFuncErr->SetMarkerSize(0);

    fAngCorrFunc->SetParameters(fMinimum->X()[0],fMinimum->X()[1],fMinimum->X()[2]);

    Int_t padpos=1;

    if(fCurrentMode[0] == '1') {
        fCanvas->cd(3)->cd(1)->cd(padpos);
        gAngCorrCorr.at(kMode_FF)->Draw("ap");
        fAngCorFuncFit.at(kMode_FF)->Draw("same");
        fAngCorFuncFit_Err.at(kMode_FF)->Draw("e3");

        padpos++;
    }
    if(fCurrentMode[1] == '1') {
        fCanvas->cd(3)->cd(1)->cd(padpos);
        gAngCorrCorr.at(kMode_FI)->Draw("ap");
        fAngCorFuncFit.at(kMode_FI)->Draw("same");
        fAngCorFuncFit_Err.at(kMode_FI)->Draw("e3");

        padpos++;
    }
    if(fCurrentMode[2] == '1') {
        fCanvas->cd(3)->cd(1)->cd(padpos);
        gAngCorrCorr.at(kMode_II)->Draw("ap");
        fAngCorFuncFit.at(kMode_II)->Draw("same");
        fAngCorFuncFit_Err.at(kMode_II)->Draw("e3");

        padpos++;
    }

    fCanvas->cd(3)->cd(2);
    fAngCorrFunc->Draw();
    fAngCorrFuncErr->Draw("e3");

    if(fTheoAngCorrFunc != nullptr) SetTransition(ftwoI1, ftwoI2, ftwoI3, fmix1, fmix2);
}

void FitQValues(bool fitQ2, bool fitQ4, double a2, double a4)
{
    //    Project();

    if(fAngCorrFunc == nullptr) {
        ERR_MESS << "Oups, you need to project you gate before"<<ENDL;
        return;
    }

    Minimize(!fitQ2,!fitQ4,true,a2,a4);

    INFO_MESS << "Use SetQValues() to load the fitted values in the program" << ENDL;
}

Double_t AngCorrFuncion(Double_t*xx,Double_t*pp) {

    Double_t x = xx[0];

    Double_t A0  = pp[0];
    Double_t A2  = pp[1];
    Double_t A4  = pp[2];

    Double_t total = A0*(1+A2*ROOT::Math::legendre(2,cos(x*TMath::DegToRad())) + A4*ROOT::Math::legendre(4,cos(x*TMath::DegToRad())));

    return total;
}

void SetTransition(int twoI1, int twoI2, int twoI3, double mix1, double mix2)
{
    if(fCanvas == nullptr) {
        ERR_MESS << "Use plot Angle and Project before"<<ENDL;
        return;
    }

    ftwoI1 = twoI1;
    ftwoI2 = twoI2;
    ftwoI3 = twoI3;
    fmix1  = mix1;
    fmix2  = mix2;

    //Determine multipolarity of gammas
    //So first find smallest possible L1 L2
    ftwoL1  = TMath::Max(abs(twoI1-twoI2),2); // no monopole
    ftwoL1p = ftwoL1+2;

    ftwoL2  = TMath::Max(abs(twoI3-twoI2),2); // no monopole
    ftwoL2p = ftwoL2+2;

    //Fix max k
    kmax = TMath::Min(ftwoI2, TMath::Min(ftwoL1p,ftwoL2p));

    if(fTheoAngCorrFunc==nullptr) {
        fTheoAngCorrFunc = new TF1("DefaultGammaCorr", Evaluate,0,180,3);
        fTheoAngCorrFunc->SetParameter(0,1);
        fTheoAngCorrFunc->SetParameter(1,0);
        fTheoAngCorrFunc->SetParameter(2,0);

        fTheoAngCorrFunc->SetLineColor(kMagenta);
    }

    if(fAngCorrFunc) fTheoAngCorrFunc->SetParameter(0,fAngCorrFunc->GetParameter(0));

    fTheoAngCorrFunc->SetParameter(1,mix1);
    fTheoAngCorrFunc->SetParameter(2,mix2);

    fCanvas->cd(3)->cd(2);
    fTheoAngCorrFunc->Draw("same");
    fTheoAngCorrFunc->Eval(0);

    Float_t _A2 = 0.;
    Float_t _A4 = 0.;

    if(Akks.size()>=1) {
        if(Akks.size()>=1) _A2 = Akks[0];
        if(Akks.size()>=2) _A4 = Akks[1];
    }
    cout << "Theoretical function: A2=" << _A2 << " A4=" << _A4 << " Delta1=" << fTheoAngCorrFunc->GetParameter(1) << ", Delta2=" << fTheoAngCorrFunc->GetParameter(2) << endl;

    Int_t padpos=1;

    for(int i=0 ; i<3 ; i++) {
        if(fCurrentMode[i] == '1') {
            fCanvas->cd(3)->cd(1)->cd(padpos);
            delete fAngCorFuncFitTheo.at(i);
            fAngCorFuncFitTheo.at(i) = new TF1(Form("AngCorFuncFit_%d",i),AngCorrFuncion,0, 180,3);
            fAngCorFuncFitTheo.at(i)->SetNpx(500);
            fAngCorFuncFitTheo.at(i)->SetLineColor(gAngCorr.at(i)->GetMarkerColor());
            fAngCorFuncFitTheo.at(i)->SetLineStyle(kDashed);
            fAngCorFuncFitTheo.at(i)->SetParameters(fMinimum->X()[0],_A2*fQValues.at(i).first.first,_A4*fQValues.at(i).second.first);
            fAngCorFuncFitTheo.at(i)->Draw("same");
            padpos++;
        }
    }
}

void FitDelta2D() {

    Float_t _A2 = 0.;
    Float_t _A4 = 0.;

    if(Akks.size()>=1) {
        if(Akks.size()>=1) _A2 = Akks[0];
        if(Akks.size()>=2) _A4 = Akks[1];
    }

    UpdateA2A42D();

    fCanvas->cd(2);

    gPad->Modified();
    gPad->Update();
}

void FitDelta(Int_t delta=1) {

    if(delta<1 || delta>2) {
        ERR_MESS << "Delta should be either 1 (first transition), or 2 (second transition)" << ENDL;
        return;
    }

    fFitMixingMode=delta;

    Float_t _A2 = 0.;
    Float_t _A4 = 0.;

    if(Akks.size()>=1) {
        if(Akks.size()>=1) _A2 = Akks[0];
        if(Akks.size()>=2) _A4 = Akks[1];
    }

    UpdateA2A4();

    if(delta==1)
        cout << "Theoretical function: A2=" << _A2 << " A4=" << _A4 << " Delta2=" << fTheoAngCorrFunc->GetParameter(2) << endl;
    else
        cout << "Theoretical function: A2=" << _A2 << " A4=" << _A4 << " Delta1=" << fTheoAngCorrFunc->GetParameter(1) << endl;


    fCanvas->cd(2);

    gPad->Modified();
    gPad->Update();

    Float_t XVal = gPad->GetFrame()->GetX2()-0.03*(gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1());
    Float_t YVal = gPad->GetFrame()->GetY2()-0.05*(gPad->GetFrame()->GetY2()-gPad->GetFrame()->GetY1());

    TString Text;
    TLatex *text;

    XVal = gPad->GetFrame()->GetX2()-0.03*(gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1());
    YVal = gPad->GetFrame()->GetY1()+0.05*(gPad->GetFrame()->GetY2()-gPad->GetFrame()->GetY1());

    Text = Form("Real #delta_{%d}: %.3f",fFitMixingMode,fTheoAngCorrFunc->GetParameter(fFitMixingMode));
    text = new TLatex(XVal,YVal,Text);
    text->SetTextColor(kGreen);text->Draw();
    text->SetTextSize(0.06);
    text->SetTextFont(132);
    text->SetTextAlign(32);
    text->Draw();

    XVal = gPad->GetFrame()->GetX1()+0.03*(gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1());
    YVal = gPad->GetFrame()->GetY1()+0.05*(gPad->GetFrame()->GetY2()-gPad->GetFrame()->GetY1());

    Text = Form("L1/L1'-L2: %d/%d-%d",ftwoL1/2,ftwoL1p/2,ftwoL2/2);
    text = new TLatex(XVal,YVal,Text);
    text->SetTextColor(kGreen);text->Draw();
    text->SetTextSize(0.06);
    text->SetTextFont(132);
    text->SetTextAlign(12);
    text->Draw();

    XVal = gPad->GetFrame()->GetX1()+0.03*(gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1());
    YVal = gPad->GetFrame()->GetY2()-0.05*(gPad->GetFrame()->GetY2()-gPad->GetFrame()->GetY1());

    TString I1=Form("%d",ftwoI1/2);
    if(ftwoI1%2==1) I1=Form("%d/2",ftwoI1);
    TString I2=Form("%d",ftwoI2/2);
    if(ftwoI2%2==1) I2=Form("%d/2",ftwoI2);
    TString I3=Form("%d",ftwoI3/2);
    if(ftwoI3%2==1) I3=Form("%d/2",ftwoI3);
    Text = Form("I1-I2-I3: %s-%s-%s",I1.Data(),I2.Data(),I3.Data());
    text = new TLatex(XVal,YVal,Text);
    text->SetTextColor(kGreen);text->Draw();
    text->SetTextSize(0.06);
    text->SetTextFont(132);
    text->SetTextAlign(12);
    text->Draw();
}

void UpdateA2A4()
{
    if(gDummy == nullptr) {
        gDummy = new TGraph;
        gDummy->SetPoint(0,-1,-1);
        gDummy->SetPoint(1,1,1);
    }

    if(gA2A4 == nullptr) {
        gA2A4 = new TGraph;
        gA2A4->GetXaxis()->SetTitle("A2");
        gA2A4->GetXaxis()->CenterTitle();
        gA2A4->GetYaxis()->SetTitle("A4");
        gA2A4->GetYaxis()->CenterTitle();
        gA2A4->SetMarkerStyle(20);
        gA2A4->SetMarkerSize(0.5);
        gA2A4->SetLineColor(kGreen);
        gA2A4->SetLineStyle(kDashed);
        gA2A4->SetLineWidth(2);
        gA2A4->GetXaxis()->SetTitleSize(0.05);
        gA2A4->GetXaxis()->SetLabelSize(0.05);
        gA2A4->GetYaxis()->SetTitleSize(0.05);
        gA2A4->GetYaxis()->SetLabelSize(0.05);
    }

    gA2A4->Set(0);

    double deltapres = fTheoAngCorrFunc->GetParameter(fFitMixingMode);

    double deltapmin=-INFINITY;
    double deltapmax=INFINITY;
    int deltanb=5000;

    double atandeltamin = atan(deltapmin);
    double atandeltamax = atan(deltapmax);

    double deltapatandelta=(atandeltamax-atandeltamin)/(deltanb+1.);
    double atandelta=atandeltamin;

    int j=0;

    Float_t _A2 = 0.;
    Float_t _A4 = 0.;

    while(atandelta<=atandeltamax) {
        fTheoAngCorrFunc->SetParameter(fFitMixingMode,tan(atandelta));
        fTheoAngCorrFunc->Eval(0);

        _A2 = 0.;
        _A4 = 0.;

        if(Akks.size()>=1) _A2 = Akks.at(0);
        if(Akks.size()>=2) _A4 = Akks.at(1);

        if(j%100==0) gA2A4->SetPoint(gA2A4->GetN(),_A2,_A4);

        atandelta += deltapatandelta;
        j++;
    }

    double minx,maxx,miny,maxy;
    minx = *min_element(gA2A4->GetX(), gA2A4->GetX()+gA2A4->GetN());
    miny = *min_element(gA2A4->GetY(), gA2A4->GetY()+gA2A4->GetN());
    maxx = *max_element(gA2A4->GetX(), gA2A4->GetX()+gA2A4->GetN());
    maxy = *max_element(gA2A4->GetY(), gA2A4->GetY()+gA2A4->GetN());

    minx = min(minx,fA2Exp-fA2Exp_errlow);
    maxx = max(maxx,fA2Exp+fA2Exp_errhigh);
    miny = min(miny,fA4Exp-fA4Exp_errlow);
    maxy = max(maxy,fA4Exp+fA4Exp_errhigh);

    double deltaX = (maxx-minx)*0.1;
    double deltaY = (maxy-miny)*0.1;

    gDummy->GetYaxis()->SetRangeUser(miny-deltaY,maxy+deltaY);
    gDummy->GetXaxis()->SetRangeUser(minx-deltaX,maxx+deltaX);

    fTheoAngCorrFunc->SetParameter(fFitMixingMode,deltapres);
    fTheoAngCorrFunc->Eval(0);

    _A2 = 0.;
    _A4 = 0.;
    if(Akks.size()>=1) _A2 = Akks.at(0);
    if(Akks.size()>=2) _A4 = Akks.at(1);

    fCanvas->cd(2);
    gDummy->Draw("ap");
    gA2A4->Draw("pc");

    if(fAngCorrFunc) {
        delete fCurrentPoint;
        fCurrentPoint = new TGraphAsymmErrors;
        fCurrentPoint->SetMarkerStyle(20);
        fCurrentPoint->SetMarkerColor(kRed);
        fCurrentPoint->SetFillStyle(3001);
        fCurrentPoint->SetFillColor(kRed);
        fCurrentPoint->SetPoint(0,fA2Exp,fA4Exp);
        fCurrentPoint->SetPointEXlow(0,fA2Exp_errlow);
        fCurrentPoint->SetPointEXhigh(0,fA2Exp_errhigh);
        fCurrentPoint->SetPointEYlow(0,fA4Exp_errlow);
        fCurrentPoint->SetPointEYhigh(0,fA4Exp_errhigh);
        fCurrentPoint->Draw("2");
        fCurrentPoint->Draw("p");

        delete fCurrentTheoPoint;
        fCurrentTheoPoint = new TGraphErrors;
        fCurrentTheoPoint->SetMarkerStyle(20);
        fCurrentTheoPoint->SetMarkerColor(kGreen);
        fCurrentTheoPoint->SetPoint(0,_A2,_A4);
        fCurrentTheoPoint->Draw("p");
    }

    // Chi2 func
    if(fAngCorrFunc) {
        delete fChi2Func;
        delete fChi2FuncDeriv;
        delete fChi2FuncDeriv2;
        fChi2Func = new TF1("Chi2Func",EvalChi2,-90,90);
        fChi2Func->SetNpx(5000);
        fChi2Func->GetXaxis()->SetTitle(Form("ATan(#delta_{%d})",fFitMixingMode));
        fChi2Func->GetYaxis()->SetTitle("#Chi^{2}/ndf");
        fChi2Func->GetXaxis()->SetTitleSize(0.05);
        fChi2Func->GetXaxis()->SetLabelSize(0.05);
        fChi2Func->GetYaxis()->SetTitleSize(0.05);
        fChi2Func->GetYaxis()->SetLabelSize(0.05);
        fChi2Func->GetXaxis()->CenterTitle();
        fChi2Func->GetYaxis()->CenterTitle();

        fChi2FuncDeriv = new TF1("Chi2FuncDeriv",EvalChi2Deriv,-90,90);
        fChi2FuncDeriv2 = new TF1("Chi2FuncDeriv2",EvalChi2Deriv2,-90,90);
        fChi2FuncDeriv->SetNpx(5000);
        fChi2FuncDeriv2->SetNpx(5000);
    }

    // Look for local minima

    double rangemin=-90;
    double rangemax=89;
    double min=0.;

    vminima.clear();

    gErrorIgnoreLevel=kFatal;
    while(true && !(_A2==0. && _A4 == 0.) && rangemin<=rangemax) {
        min = fChi2FuncDeriv->GetX(0.,rangemin,rangemax);
        if(isnan(min)) break;
        if(fChi2FuncDeriv2->Eval(min)>0) {
            vminima.push_back(min);
        }
        rangemin = min+0.5;
    }
    gErrorIgnoreLevel=kInfo;

    fCanvas->cd(4);
    gPad->SetLogy(1);
    gPad->SetLogz(0);

    fChi2Func->Draw();

    gPad->Modified();
    gPad->Update();

    double BestChi2 = numeric_limits<double>::max();
    BestDelta=0.;

    int nmin=0;
    for(auto min: vminima) {
        fCanvas->cd(4);
        gPad->Modified();
        gPad->Update();

        double y1 = TMath::Power(10,gPad->GetUymin());
        double y2 = TMath::Power(10,gPad->GetUymax());

        Color_t col = kBlue;
        if(nmin==1) col = kMagenta;
        if(nmin==2) col = kCyan;
        nmin++;

        TLine *l = new TLine(min,y1,min,y2);
        l->SetLineColor(col);
        l->SetLineStyle(kDashed);
        l->Draw();

        TLatex *text = new TLatex(min-0.005*(gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1()),y2,Form("#delta_{%d}=%.3f, #Chi^{2}=%.1f ",fFitMixingMode,tan(min*TMath::DegToRad()),fChi2Func->Eval(min)));
        text->SetTextColor(col);
        text->SetTextSize(0.05);
        text->SetTextFont(132);
        text->SetTextAlign(31);
        text->SetTextAngle(90);
        text->Draw();

        fCanvas->cd(2);

        double deltapres = fTheoAngCorrFunc->GetParameter(fFitMixingMode);
        fTheoAngCorrFunc->SetParameter(fFitMixingMode,tan(min*TMath::DegToRad()));
        fTheoAngCorrFunc->Eval(0);

        _A2 = 0.;
        _A4 = 0.;
        if(Akks.size()>=1) _A2 = Akks.at(0);
        if(Akks.size()>=2) _A4 = Akks.at(1);

        TGraph *gr = new TGraph;
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(col);

        gr->SetPoint(0,_A2,_A4);
        gr->Draw("p");

        if(fChi2Func->Eval(min)<BestChi2) {
            BestChi2=fChi2Func->Eval(min);
            BestDelta = tan(min*TMath::DegToRad());
        }

        fTheoAngCorrFunc->SetParameter(fFitMixingMode,deltapres);
        fTheoAngCorrFunc->Eval(0);
    }
}

void UpdateA2A42D()
{
    if(gDummy == nullptr) {
        gDummy = new TGraph;
        gDummy->SetPoint(0,-1,-1);
        gDummy->SetPoint(1,1,1);
    }

    if(gA2A4 == nullptr) {
        gA2A4 = new TGraph;
        gA2A4->GetXaxis()->SetTitle("A2");
        gA2A4->GetXaxis()->CenterTitle();
        gA2A4->GetYaxis()->SetTitle("A4");
        gA2A4->GetYaxis()->CenterTitle();
        gA2A4->SetMarkerStyle(1);
        gA2A4->SetMarkerSize(0.5);
        gA2A4->SetLineColor(kGreen);
        gA2A4->SetLineStyle(kDashed);
        gA2A4->SetLineWidth(2);
        gA2A4->GetXaxis()->SetTitleSize(0.05);
        gA2A4->GetXaxis()->SetLabelSize(0.05);
        gA2A4->GetYaxis()->SetTitleSize(0.05);
        gA2A4->GetYaxis()->SetLabelSize(0.05);
    }

    gA2A4->Set(0);

    double deltapres1 = fTheoAngCorrFunc->GetParameter(1);
    double deltapres2 = fTheoAngCorrFunc->GetParameter(2);

    double deltapmin=-INFINITY;
    double deltapmax=INFINITY;
    int deltanb=5000;

    double atandeltamin = atan(deltapmin);
    double atandeltamax = atan(deltapmax);

    double deltapatandelta=(atandeltamax-atandeltamin)/(deltanb+1.);
    double atandelta1=atandeltamin;

    int j=0;

    Float_t _A2 = 0.;
    Float_t _A4 = 0.;

    while(atandelta1<=atandeltamax) {
        double atandelta2=atandeltamin;
        while(atandelta2<=atandeltamax) {
            fTheoAngCorrFunc->SetParameter(1,tan(atandelta1));
            fTheoAngCorrFunc->SetParameter(2,tan(atandelta2));
            fTheoAngCorrFunc->Eval(0);

            _A2 = 0.;
            _A4 = 0.;

            if(Akks.size()>=1) _A2 = Akks.at(0);
            if(Akks.size()>=2) _A4 = Akks.at(1);

            if(j%200==0) gA2A4->SetPoint(gA2A4->GetN(),_A2,_A4);
            j++;
            atandelta2 += deltapatandelta;
        }
        atandelta1 += deltapatandelta;
    }

    double minx,maxx,miny,maxy;
    minx = *min_element(gA2A4->GetX(), gA2A4->GetX()+gA2A4->GetN());
    miny = *min_element(gA2A4->GetY(), gA2A4->GetY()+gA2A4->GetN());
    maxx = *max_element(gA2A4->GetX(), gA2A4->GetX()+gA2A4->GetN());
    maxy = *max_element(gA2A4->GetY(), gA2A4->GetY()+gA2A4->GetN());

    minx = min(minx,fA2Exp-fA2Exp_errlow);
    maxx = max(maxx,fA2Exp+fA2Exp_errhigh);
    miny = min(miny,fA4Exp-fA4Exp_errlow);
    maxy = max(maxy,fA4Exp+fA4Exp_errhigh);

    double deltaX = (maxx-minx)*0.1;
    double deltaY = (maxy-miny)*0.1;

    gDummy->GetYaxis()->SetRangeUser(miny-deltaY,maxy+deltaY);
    gDummy->GetXaxis()->SetRangeUser(minx-deltaX,maxx+deltaX);

    fTheoAngCorrFunc->SetParameter(1,deltapres1);
    fTheoAngCorrFunc->SetParameter(2,deltapres2);

    fTheoAngCorrFunc->Eval(0);

    _A2 = 0.;
    _A4 = 0.;
    if(Akks.size()>=1) _A2 = Akks.at(0);
    if(Akks.size()>=2) _A4 = Akks.at(1);

    fCanvas->cd(2);
    gDummy->Draw("ap");
    gA2A4->Draw("p");

    if(fAngCorrFunc) {
        delete fCurrentPoint;
        fCurrentPoint = new TGraphAsymmErrors;
        fCurrentPoint->SetMarkerStyle(20);
        fCurrentPoint->SetMarkerColor(kRed);
        fCurrentPoint->SetFillStyle(3001);
        fCurrentPoint->SetFillColor(kRed);
        fCurrentPoint->SetPoint(0,fA2Exp,fA4Exp);
        fCurrentPoint->SetPointEXlow(0,fA2Exp_errlow);
        fCurrentPoint->SetPointEXhigh(0,fA2Exp_errhigh);
        fCurrentPoint->SetPointEYlow(0,fA4Exp_errlow);
        fCurrentPoint->SetPointEYhigh(0,fA4Exp_errhigh);
        fCurrentPoint->Draw("2");
        fCurrentPoint->Draw("p");

        delete fCurrentTheoPoint;
        fCurrentTheoPoint = new TGraphErrors;
        fCurrentTheoPoint->SetMarkerStyle(20);
        fCurrentTheoPoint->SetMarkerColor(kGreen);
        fCurrentTheoPoint->SetPoint(0,_A2,_A4);
        fCurrentTheoPoint->Draw("p");
    }

    // Chi2 func
    if(fAngCorrFunc) {
        delete fChi2Func2D;
//        delete fChi2FuncDeriv;
//        delete fChi2FuncDeriv2;
        fChi2Func2D = new TH2F("Chi2Func","Chi2Func",180,-90,90,180,-90,90);
        fChi2Func2D->GetXaxis()->SetTitle(Form("ATan(#delta_{%d})",1));
        fChi2Func2D->GetYaxis()->SetTitle(Form("ATan(#delta_{%d})",2));
        fChi2Func2D->GetZaxis()->SetTitle("#Chi^{2}/ndf");
        fChi2Func2D->GetXaxis()->SetTitleSize(0.05);
        fChi2Func2D->GetXaxis()->SetLabelSize(0.05);
        fChi2Func2D->GetYaxis()->SetTitleSize(0.05);
        fChi2Func2D->GetYaxis()->SetLabelSize(0.05);
        fChi2Func2D->GetXaxis()->CenterTitle();
        fChi2Func2D->GetYaxis()->CenterTitle();

//        fChi2FuncDeriv = new TF1("Chi2FuncDeriv",EvalChi2Deriv,-90,90);
//        fChi2FuncDeriv2 = new TF1("Chi2FuncDeriv2",EvalChi2Deriv2,-90,90);
//        fChi2FuncDeriv->SetNpx(5000);
//        fChi2FuncDeriv2->SetNpx(5000);
    }

    double BestChi2 = numeric_limits<double>::max();
    double BestDelta1=0;
    double BestDelta2=0;
    double BestTanDelta1=0;
    double BestTanDelta2=0;

    for(int i=1 ; i<fChi2Func2D->GetNbinsX() ; i++) {
        for(int j=1 ; j<fChi2Func2D->GetNbinsY() ; j++) {
            double atan_delta1 = fChi2Func2D->GetXaxis()->GetBinLowEdge(i);
            double atan_delta2 = fChi2Func2D->GetYaxis()->GetBinLowEdge(j);
            double Chi2 = EvalChi2_2D(atan_delta1,atan_delta2);
            if(Chi2<=BestChi2) {
                BestChi2 = Chi2;
                BestTanDelta1 = atan_delta1;
                BestTanDelta2 = atan_delta2;
                BestDelta1 = tan(atan_delta1*TMath::DegToRad());
                BestDelta2 = tan(atan_delta2*TMath::DegToRad());
            }
            fChi2Func2D->SetBinContent(i,j,Chi2);
        }
    }

    fCanvas->cd(4);
    gPad->SetLogy(0);
    fChi2Func2D->Draw("colz");


    TGraph *gr = new TGraph;
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kGreen);
    gr->SetPoint(0,BestTanDelta1,BestTanDelta2);
    gr->Draw("p");

    TLatex *text = new TLatex(BestTanDelta1,BestTanDelta2-5,Form("#delta_{1}=%.3f, #delta_{2}=%.3f, #Chi^{2}=%.1f ",BestDelta1,BestDelta2,BestChi2));
    text->SetTextColor(kGreen);
    text->SetTextSize(0.05);
    text->SetTextFont(132);
    text->SetTextAlign(23);
    text->SetTextAngle(0);
    text->Draw();

    gPad->Modified();
    gPad->Update();

    fCanvas->cd(2);

    deltapres1 = fTheoAngCorrFunc->GetParameter(1);
    deltapres2 = fTheoAngCorrFunc->GetParameter(2);

    fTheoAngCorrFunc->SetParameter(1,BestDelta1);
    fTheoAngCorrFunc->SetParameter(2,BestDelta1);
    fTheoAngCorrFunc->Eval(0);

    _A2 = 0.;
    _A4 = 0.;
    if(Akks.size()>=1) _A2 = Akks.at(0);
    if(Akks.size()>=2) _A4 = Akks.at(1);

    gr = new TGraph;
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);

    gr->SetPoint(0,_A2,_A4);
    gr->Draw("p");

    fTheoAngCorrFunc->SetParameter(1,deltapres1);
    fTheoAngCorrFunc->SetParameter(2,deltapres2);
    fTheoAngCorrFunc->Eval(0);

    gPad->Modified();
    gPad->Update();
}


// Thesis Marcelo Barbosa, p 40
// https://cds.cern.ch/record/1641449/files/CERN-THESIS-2010-294.pdf?subformat=pdfa
// Fk(LL′IiI)
double Fk(int twoL1, int twoL1p, int twoIi, int twoI, int k)
{
    double tot = pow(-1.,(twoIi+twoI-2)/2)*
                 sqrt((twoL1+1)*(twoL1p+1)*(twoI+1)*(2*k+1))*
                 ROOT::Math::wigner_3j(twoL1,twoL1p,2*k,2,-2,0)*          // wigner_3j takes 2ji as input
                 ROOT::Math::wigner_6j(twoL1,twoL1p,2*k,twoI,twoI,twoIi); // wigner_6j takes 2ji as input

    return  tot;
}

void Eval_Ak()
{
    double Delta1 = fTheoAngCorrFunc->GetParameter(1);
    double Delta2 = fTheoAngCorrFunc->GetParameter(2);

    Akks.clear();

    for(int k=2; k<=kmax; k+=2) {
        double Ak1 = 1./(1.+Delta1*Delta1)*(Fk(ftwoL1,ftwoL1,ftwoI1,ftwoI2,k) - 2.*Delta1*Fk(ftwoL1,ftwoL1p,ftwoI1,ftwoI2,k) + Delta1*Delta1*Fk(ftwoL1p,ftwoL1p,ftwoI1,ftwoI2,k));
        double Ak2 = 1./(1.+Delta2*Delta2)*(Fk(ftwoL2,ftwoL2,ftwoI3,ftwoI2,k) + 2.*Delta2*Fk(ftwoL2,ftwoL2p,ftwoI3,ftwoI2,k) + Delta2*Delta2*Fk(ftwoL2p,ftwoL2p,ftwoI3,ftwoI2,k));
        double Akk = Ak1*Ak2;
        Akks.push_back(Akk);
    }
}

double Evaluate(double *x,double *parameters)
{
    // W(theta) = 1 + Sum(Ak*P(cos(thehta)))

    Eval_Ak();

    double tot=1;
    int iak=0;
    for(int k=2; k<=kmax; k+=2, iak++) {
        double arg = cos(*x*TMath::DegToRad());
        double qcor=1.;
        //        if(k==2) qcor=Q2;
        //        if(k==4) qcor=Q4;
        tot+=Akks.at(iak)*qcor*ROOT::Math::legendre(k,arg);
    }
    return parameters[0]*tot;
}

double EvalChi2(double *x,double */*parameters*/)
{
    // From Urban_2013_J_Inst_8_P03014

    if(fTheoAngCorrFunc==nullptr || fAngCorrFunc == nullptr) return 0.;

    double deltapres = fTheoAngCorrFunc->GetParameter(fFitMixingMode);
    double delta = tan(x[0]*TMath::DegToRad());

    fTheoAngCorrFunc->SetParameter(fFitMixingMode,delta);
    fTheoAngCorrFunc->Eval(0);

    double Chi2_A2=0;
    double Chi2_A4=0;

    double _A2 = 0.;
    double _A4 = 0.;

    if(Akks.size()>=1) _A2 = Akks.at(0);
    if(Akks.size()>=2) _A4 = Akks.at(1);

    Chi2_A2 = (fAngCorrFunc->GetParameter(1)-_A2)/fAngCorrFunc->GetParError(1);
    Chi2_A4 = (fAngCorrFunc->GetParameter(2)-_A4)/fAngCorrFunc->GetParError(2);

    fTheoAngCorrFunc->SetParameter(fFitMixingMode,deltapres);
    fTheoAngCorrFunc->Eval(0);

    return Chi2_A2*Chi2_A2+Chi2_A4*Chi2_A4;
}

double EvalChi2_2D(double atandelta1,double atandelta2)
{
    // From Urban_2013_J_Inst_8_P03014

    if(fTheoAngCorrFunc==nullptr || fAngCorrFunc == nullptr) return 0.;

    double deltapres1 = fTheoAngCorrFunc->GetParameter(1);
    double deltapres2 = fTheoAngCorrFunc->GetParameter(2);

    double delta1 = tan(atandelta1*TMath::DegToRad());
    double delta2 = tan(atandelta2*TMath::DegToRad());

    fTheoAngCorrFunc->SetParameter(1,delta1);
    fTheoAngCorrFunc->SetParameter(2,delta2);
    fTheoAngCorrFunc->Eval(0);

    double Chi2_A2=0;
    double Chi2_A4=0;

    double _A2 = 0.;
    double _A4 = 0.;

    if(Akks.size()>=1) _A2 = Akks.at(0);
    if(Akks.size()>=2) _A4 = Akks.at(1);

    Chi2_A2 = (fAngCorrFunc->GetParameter(1)-_A2)/fAngCorrFunc->GetParError(1);
    Chi2_A4 = (fAngCorrFunc->GetParameter(2)-_A4)/fAngCorrFunc->GetParError(2);

    fTheoAngCorrFunc->SetParameter(1,deltapres1);
    fTheoAngCorrFunc->SetParameter(2,deltapres2);
    fTheoAngCorrFunc->Eval(0);

    return Chi2_A2*Chi2_A2+Chi2_A4*Chi2_A4;
}

double EvalChi2Deriv(double *x,double */*parameters*/)
{
    return fChi2Func->Derivative(x[0]);
}

double EvalChi2Deriv2(double *x,double */*parameters*/)
{
    return fChi2Func->Derivative2(x[0]);
}

void SaveCanvas() {
    if(fCanvas == nullptr) return;

    TString Name = Form("AngCorr_%.1f_%.1f.png",GatePx.first,GatePy.first);
    fCanvas->SaveAs(Name);
}
