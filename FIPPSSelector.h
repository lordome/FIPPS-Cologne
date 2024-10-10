//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  9 17:18:32 2019 by ROOT version 6.18/04
// from TChain FIPPS_Tree/
//////////////////////////////////////////////////////////

#ifndef FIPPSSelector_h
#define FIPPSSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TEnv.h>
#include <TParameter.h>
#include <TH2.h>
#include <TH3.h>
#include <TVector3.h>

#include <progressbar.h>
#include <array>
#include <tuple>

using namespace std;

class FIPPSSelector : public TSelector {

    struct GGEffAngCorStruc {
        Float_t  EGate;
        Float_t  EWidth;
        Float_t  a2;
        Float_t  a4;
    };

    struct GGIsomerStruct {
        Float_t  G1Min;
        Float_t  G1Max;
        Float_t  G2Min;
        Float_t  G2Max;

        TH2 *hist = nullptr;

        Bool_t IsValid(Long64_t Tgamma1, Long64_t Tgamma2) {
            return( ((Tgamma1>=G1Min && Tgamma1<=G1Max) && (Tgamma2>=G2Min && Tgamma2<=G2Max)) ||
                    ((Tgamma2>=G1Min && Tgamma2<=G1Max) && (Tgamma1>=G2Min && Tgamma1<=G2Max))
                    );
        }
    };

public:

    TTree *fCurrentTree = nullptr;

    Bool_t fUseProof    = false;
    Int_t fMaxADC       = 0;
    Int_t fMaxIds       = 0;
    Int_t fMaxNClov     = 0;
    Int_t fNGeDet       = 0;
    Int_t fNFIPPS_GeDet = 0;
    Int_t fNIFIN_GeDet  = 0;
    Int_t fNLaBr3	= 0;
    Int_t fNTAG         = 0;
    Int_t fNAC          = 0;

    Int_t fCurrentRun   = 0;
    Int_t fReferenceRun = 0;
    TDatime fFirstDate;
    TString fWorkingDir;

    Int_t fRunNr;
    Bool_t fUseTAG = false;
    Bool_t fUseLaBr3 = false;

    Bool_t fPlotBasic = false;
    Bool_t fPlotMonitoring = false;
    Bool_t fPlotRates = false;
    Bool_t fPlotCalibrated = false;
    Bool_t fPlotRelAngles = false;
    Bool_t fPlotXTalk = false;
    Bool_t fPlotDevHists = false;
    Int_t fEWindowEnergy = 0;
    Int_t fEWindowWidth = 0;

    // to plot TSC matrix
    Bool_t fPlotOsloTSC  = false;
    
    // to set an energy threshold on coincidences for GG (ex. E1 + E2 > 3MeV)
    Bool_t fUseSumThresholdForGG = false;
    TString fSumThresholdForGG ;

    
    TObjArray *fArrayOfPositions = nullptr;

    //if not using proof
    TFile *fFinalFile;
    TString fFileNameOut;

    Float_t fRelAngleBin_Width;

    vector< pair<Float_t, Float_t > > fRelAnglesGates;

    // Gate/Theta/Mode
    vector< map < Int_t, array < TH2*, 3 > > > fRelativeAnglesGG;
    vector< map < Int_t, array < TH2*, 3 > > > fRelativeAnglesGG_Norm;

/*    vector< map < Int_t, map< Int_t, TH2*> > > fRelativeAnglesGG;
    vector< map < Int_t, map< Int_t, TH2*> > > fRelativeAnglesGG_Norm*/;

    TH2* fDummyMatrix = nullptr;

    vector < pair < Int_t, Int_t > > fRelAnglesNFired; // To know if each gate has been fired and how many times: NFired, LastClovId (if NFired==1, LastClovId must not be filled)

    Int_t fDTClover;
    Int_t fTimeWindow;

    Float_t fRateScaleFactor;

    progressbar *progress;

    // Histograms
    TH1 *hGlob_Mult = nullptr;
    TH1 *hMeanMultByEvt = nullptr;
    TH1 *hRelativeAnglesDist = nullptr;

    TH1 *hADC_Values = nullptr;
    TH1 *hIds = nullptr;
    TH2 *hADC_vs_Ids = nullptr;
    TH2 *hIds_vs_Ids = nullptr;

    TH2 *hTiming_Window_DT = nullptr;
    TH2 *hTiming_DT_to_Id0 = nullptr;
    TH2 *hTiming_DT_vsEnergyTAG = nullptr;
    TH2 *hTiming_Clov_TAG = nullptr;


    TH1 *hGe_ESpectra[200];
    TH1 *hGe_ESpectraMult[2][200];
    TH2 *hGe_EGGM2[200][4][4];

    TH1 *hLaBr3_ESpectra[200];

    TH1 *hGe_E = nullptr;
    TH1 *hGe_E_FIPPS = nullptr;
    TH1 *hGe_E_IFIN = nullptr;

    TH2 *hGe_EvsId = nullptr;
    TH1 *hSumE[200];

    TH2 *hTAG_vs_TAG = nullptr;

    TH1 *hAddE[200];
    TH1 *hERej[200];

    TH2 *hEGammavsId = nullptr;
    TH1 *hEGamma = nullptr;
    TH1 *hEGamma_FIPPS = nullptr;
    TH1 *hEGamma_IFIN = nullptr;
    TH1 *hEGamma_FTag = nullptr;
    TH1 *hEGamma_BTag = nullptr;
    TH2 *hEGamma_Qtot = nullptr;
    TH2 *hEGamma_FTag_QTot = nullptr;
    TH2 *hDeltaT_EGamma =nullptr;
    TH2 *hEGG_FTag = nullptr;
    TH2 *hEGG_BTag = nullptr;
    TH2 *hEGG = nullptr;
    TH2 *hEGG_ESumThreshold = nullptr;

    TH2 *hEGamma_ETot = nullptr;
    TH2 *hEGamma_ETot_M2 = nullptr;
    TH2 *hGe_E_GG = nullptr;
    TH2 *hGe_E_GG_FIPPS_FIPPS = nullptr;
    TH2 *hGe_E_GG_IFIN_IFIN = nullptr;
    TH2 *hGe_E_GG_FIPPS_IFIN = nullptr;
    TH2 *hGe_E_GG_FIPPS_IFIN_sym = nullptr;

    TH2 *hEGammavsTime = nullptr;
    TH2 *hEGammaCalvsRun[200];
    TH2 *hQTotvsTime[200];
    TH2 *hQTotvsRun[200];
    TH1 *histoQTot[200];
    TH1 *hFissionRate = nullptr;
    TH1 *hACRate = nullptr;
    TH1 *hFIPPSRate = nullptr;
    TH1 *hIFINRate = nullptr;
    TH1 *hTAGRate = nullptr;
    TH1 *hADC45Rate = nullptr;

    map <Int_t, map< Int_t, TH2*>> hMapOfClovbyClovGG;

    bool fDoAngularUsingCloverAngles = false;
    bool fDoAngularCorrelations = false;
    bool fDoNormalizarion = false;

    pair< tuple<int, float, float>, tuple<int, float, float> > fEGammaETotBinning;
    pair< tuple<int, float, float>, tuple<int, float, float> > fRelAnglesBinning;

    tuple<int, float, float> fEGBinning;
    tuple<int, float, float> fEGGBinning;


    static const int NLastVecs = 5;
    Int_t NFilled=0;

    TString fGGIsomers;

    TString fGGEff_AngCorrCond;
    vector < GGEffAngCorStruc > fListofGGEffAngCorrCond;
    vector < array <TH1*,3> > fListofGGEffAngCorrCondHists;
    vector < array <TH1*,3> > fListofGGAddEffAngCorrCondHists;

    vector < GGIsomerStruct > fListofGGIsomers;

    // keep last NLastVecs positions (for X and Y binning)
    vector< vector< tuple<Float_t, Float_t, Float_t, Int_t> > > fLastPos_BinX;
    vector< vector< tuple<Float_t, Float_t, Float_t, Int_t> > > fLastPos_BinY;

    array< pair<TVector3, TVector3>, NLastVecs > fLastPositions;

public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
    TTreeReaderValue<Int_t> multiplicity = {fReader, "multiplicity"};
    TTreeReaderArray<Int_t> energy = {fReader, "energy"};
    TTreeReaderArray<Int_t> extraData = {fReader, "extraData"};
    TTreeReaderArray<Short_t> adc = {fReader, "adc"};
    TTreeReaderArray<Long64_t> timeStamp = {fReader, "timeStamp"};
    TTreeReaderArray<Short_t> globid = {fReader, "globid"};
    TTreeReaderArray<Short_t> type = {fReader, "type"};

    TTreeReaderValue<Int_t> Ge_Mult = {fReader, "Ge_Mult"};
    TTreeReaderArray<Short_t> Ge_Id = {fReader, "Ge_Id"};
    TTreeReaderArray<UChar_t> Ge_Type = {fReader, "Ge_Type"};
    TTreeReaderArray<UChar_t> Ge_ClovMult = {fReader, "Ge_ClovMult"};
    TTreeReaderArray<Short_t> Ge_ClovId = {fReader, "Ge_ClovId"};
    TTreeReaderArray<Float_t> Ge_E = {fReader, "Ge_E"};
    TTreeReaderArray<Float_t> Ge_Extra = {fReader, "Ge_Extra"};
    TTreeReaderArray<Float_t> Ge_X = {fReader, "Ge_X"};
    TTreeReaderArray<Float_t> Ge_Y = {fReader, "Ge_Y"};
    TTreeReaderArray<Float_t> Ge_Z = {fReader, "Ge_Z"};
    TTreeReaderArray<Long64_t> Ge_TS = {fReader, "Ge_TS"};

    TTreeReaderValue<Int_t> Add_Mult = {fReader, "Add_Mult"};
    TTreeReaderArray<Short_t> Add_ClovId = {fReader, "Add_ClovId"};
    TTreeReaderArray<Float_t> Add_E = {fReader, "Add_E"};
    TTreeReaderArray<Float_t> Add_ERej = {fReader, "Add_ERej"};
    TTreeReaderArray<Float_t> Add_X = {fReader, "Add_X"};
    TTreeReaderArray<Float_t> Add_Y = {fReader, "Add_Y"};
    TTreeReaderArray<Float_t> Add_Z = {fReader, "Add_Z"};
    TTreeReaderArray<Long64_t> Add_TS = {fReader, "Add_TS"};

    Int_t fLast_AddMult[NLastVecs];
    Float_t fLast_Add_E[NLastVecs][100];
    Float_t fLast_Add_X[NLastVecs][100];
    Float_t fLast_Add_Y[NLastVecs][100];
    Float_t fLast_Add_Z[NLastVecs][100];
    Long64_t fLast_Add_TS[NLastVecs][100];
    Short_t fLast_Add_ClovId[NLastVecs][100];

    TTreeReaderValue<Int_t> AC_Mult = {fReader, "AC_Mult"};
    TTreeReaderArray<Short_t> AC_Id = {fReader, "AC_Id"};
    TTreeReaderArray<UChar_t> AC_Type = {fReader, "AC_Type"};
    TTreeReaderArray<Short_t> AC_ClovId = {fReader, "AC_ClovId"};
    TTreeReaderArray<Float_t> AC_E = {fReader, "AC_E"};
    TTreeReaderArray<Long64_t> AC_TS = {fReader, "AC_TS"};

    /*
   TTreeReaderValue<Int_t> LaBr3_Mult = {fReader, "LaBr3_Mult"};
   TTreeReaderArray<Short_t> LaBr3_Id = {fReader, "LaBr3_Id"};
   TTreeReaderArray<Float_t> LaBr3_E = {fReader, "LaBr3_E"};
   TTreeReaderArray<Long64_t> LaBr3_TS = {fReader, "LaBr3_TS"};
*/
    std::unique_ptr<TTreeReaderValue<Int_t>> LaBr3_Mult_ptr;
    std::unique_ptr<TTreeReaderArray<Short_t>> LaBr3_Id_ptr;
    std::unique_ptr<TTreeReaderArray<Float_t>> LaBr3_E_ptr;
    std::unique_ptr<TTreeReaderArray<Long64_t>> LaBr3_TS_ptr;

    Int_t    *LaBr3_Mult;
    Float_t  *LaBr3_E;
    Short_t  *LaBr3_Id;
    Long64_t *LaBr3_TS;

    std::unique_ptr<TTreeReaderValue<Int_t>> TAG_Mult_ptr;
    std::unique_ptr<TTreeReaderArray<Float_t>> TAG_QTot_ptr;
    std::unique_ptr<TTreeReaderArray<Float_t>> TAG_QShort_ptr;
    std::unique_ptr<TTreeReaderArray<Bool_t>> TAG_IsFission_ptr;
    std::unique_ptr<TTreeReaderArray<Short_t>> TAG_Id_ptr;
    std::unique_ptr<TTreeReaderArray<Long64_t>> TAG_TS_ptr;

    Int_t    *TAG_Mult;
    Float_t  *TAG_QTot;
    Float_t  *TAG_QShort;
    Bool_t   *TAG_IsFission;
    Short_t  *TAG_Id;
    Long64_t *TAG_TS;

    TTreeReaderValue<Double_t> GlobTime = {fReader, "GlobTime"};

    FIPPSSelector(TTree * /*tree*/ =0) { }
    virtual ~FIPPSSelector() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();

    TH1 *BuildHist1D(TString HistName, Int_t NBins, Float_t Min, Float_t Max, TString AxisName, TString FolderName, TString Type="TH1F", Bool_t Timing = false);
    TH2 *BuildHist2D(TString HistName, Int_t NBinsX, Float_t MinX, Float_t MaxX, TString AxisNameX, Int_t NBinsY, Float_t MinY, Float_t MaxY, TString AxisNameY, TString FolderName, Bool_t Timing=false);
    TH3 *BuildHist3D(TString HistName, Int_t NBinsX, Float_t MinX, Float_t MaxX, TString AxisNameX, Int_t NBinsY, Float_t MinY, Float_t MaxY, TString AxisNameY, Int_t NBinsZ, Float_t MinZ, Float_t MaxZ, TString AxisNameZ, TString FolderName);


    TString GetEnvValue(TString name, TString defaulvalue){
        if(fUseProof) return ((TString)((TParameter<TString>*)fInput->FindObject(name))->GetVal());
        else return gEnv->GetValue(name,defaulvalue);
    }

    Int_t GetEnvValue(TString name, Int_t defaulvalue) {
        if(fUseProof) {
            if(fInput->FindObject(name))
                return ((Int_t)((TParameter<Int_t>*)fInput->FindObject(name))->GetVal());
            else
                return defaulvalue;
        }
        else return gEnv->GetValue(name,defaulvalue);
    }

    Double_t GetEnvValue(TString name, Double_t defaulvalue){
        if(fUseProof) return ((Double_t)((TParameter<Double_t>*)fInput->FindObject(name))->GetVal());
        else return gEnv->GetValue(name,defaulvalue);
    }

    Long64_t GetEnvValue(TString name, Long64_t defaulvalue){
        if(fUseProof) return ((Long64_t)((TParameter<Long64_t>*)fInput->FindObject(name))->GetVal());
        else return gEnv->GetValue(name,(Double_t)defaulvalue);
    }

    void SetPositions(TObjArray *arr) {fArrayOfPositions = (TObjArray*)arr->Clone();}

    ClassDef(FIPPSSelector,0);

};
#endif

#ifdef FIPPSSelector_cxx
void FIPPSSelector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    fReader.SetTree(tree);
}

Bool_t FIPPSSelector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef FIPPSSelector_cxx
