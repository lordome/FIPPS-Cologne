#define FIPPSSelector_cxx
#include "FIPPSSelector.h"

#include "TH2F.h"
#include "TH1F.h"

#include "TStyle.h"
#include "FIPPS_Softs.h"
#include "TSystem.h"
#include "TProofServ.h"
#include "TVector3.h"

void FIPPSSelector::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    if(fInput && fInput->FindObject("UseProof"))
        fUseProof = true;
    else
        fUseProof = false;

    if(fUseProof == false) {
        fFileNameOut = GetEnvValue("FileNameOut",(TString)0);
        fFinalFile = new TFile(fFileNameOut,"recreate");

        progress = new progressbar(GetEnvValue("NEntries",(Long64_t)0));
    }

    TString mode = GetEnvValue("Mode",(TString)0);
}

void FIPPSSelector::SlaveBegin(TTree * /*tree*/)
{
    TString option = GetOption();

    if(fInput && fInput->FindObject("UseProof"))
        fUseProof = true;
    else
        fUseProof = false;


    fMaxADC        = GetEnvValue("NMaxADC",(Int_t)0)+1;
    fMaxIds        = GetEnvValue("NMaxIds",(Int_t)0)+1;
    fMaxNClov      = GetEnvValue("NMaxClov",(Int_t)0)+1;
    fNFIPPS_GeDet  = GetEnvValue("NFIPPS_Ge",(Int_t)0);
    fNIFIN_GeDet   = GetEnvValue("NIFIN_Ge",(Int_t)0);
    fNGeDet        = fNFIPPS_GeDet + fNIFIN_GeDet;
    fNTAG          = GetEnvValue("NTAG",(Int_t)0);
    fNAC           = fNFIPPS_GeDet/4*3 + fNIFIN_GeDet/4;
    fNLaBr3   	   = GetEnvValue("NLaBr3",(Int_t)0);

    fDTClover      = GetEnvValue("DTClover",(Int_t)0);
    fTimeWindow    = GetEnvValue("TimeWindow",(Int_t)0);


    TString PlotOptions = GetEnvValue("PlotOptions",(TString)0);
    TString RelAnglesMode = GetEnvValue("RelAnglesMode",(TString)0);
    TString RelAnglesNorm = GetEnvValue("RelAnglesNorm",(TString)0);
    TString RelAnglesCond = GetEnvValue("RelAnglesCond",(TString)0);
    TString RelAnglesBinning = GetEnvValue("RelAnglesBinning",(TString)0);

    TString ETotBinning = GetEnvValue("EGammaETotBinning",(TString)0);

    TObjArray *arr = ETotBinning.Tokenize(" ");
    if(arr->GetEntries() == 3) {
        get<0>(fEGammaETotBinning.first) = ((TString)arr->At(0)->GetName()).Atoi();
        get<1>(fEGammaETotBinning.first) = ((TString)arr->At(1)->GetName()).Atof();
        get<2>(fEGammaETotBinning.first) = ((TString)arr->At(2)->GetName()).Atof();

        get<0>(fEGammaETotBinning.second) = ((TString)arr->At(0)->GetName()).Atoi();
        get<1>(fEGammaETotBinning.second) = ((TString)arr->At(1)->GetName()).Atof();
        get<2>(fEGammaETotBinning.second) = ((TString)arr->At(2)->GetName()).Atof();
    }
    else if(arr->GetEntries() == 6) {
        get<0>(fEGammaETotBinning.first) = ((TString)arr->At(0)->GetName()).Atoi();
        get<1>(fEGammaETotBinning.first) = ((TString)arr->At(1)->GetName()).Atof();
        get<2>(fEGammaETotBinning.first) = ((TString)arr->At(2)->GetName()).Atof();
        get<0>(fEGammaETotBinning.second) = ((TString)arr->At(3)->GetName()).Atoi();
        get<1>(fEGammaETotBinning.second) = ((TString)arr->At(4)->GetName()).Atof();
        get<2>(fEGammaETotBinning.second) = ((TString)arr->At(5)->GetName()).Atof();
    }
    else {
        ERR_MESS << "Bad EGammaETot binning definition, needs 3 (symetric) or 6 (asymetric) parameters, 4096 0 4094 used" << ENDL;
        get<0>(fEGammaETotBinning.first)  = 4096;
        get<1>(fEGammaETotBinning.first)  = 0;
        get<2>(fEGammaETotBinning.first)  = 4096;
        get<0>(fEGammaETotBinning.second) = 4096;
        get<1>(fEGammaETotBinning.second) = 0;
        get<2>(fEGammaETotBinning.second) = 4096;
    }
    delete arr;

    TString EGBinning = GetEnvValue("EGBinning",(TString)0);
    TObjArray *arrEG = EGBinning.Tokenize(" ");

    INFO_MESS << 	arrEG->GetEntries() << "Entries" << ENDL;
    if(arrEG->GetEntries() == 3) {
        get<0>(fEGBinning) = ((TString)arrEG->At(0)->GetName()).Atoi();
        get<1>(fEGBinning) = ((TString)arrEG->At(1)->GetName()).Atof();
        get<2>(fEGBinning) = ((TString)arrEG->At(2)->GetName()).Atof();
    }
    else {
        ERR_MESS << "Bad EGBinning binning definition, needs 3 parameters, 8192 0 8192 used" << ENDL;
        get<0>(fEGGBinning)  = 11000;
        get<1>(fEGGBinning)  = 0;
        get<2>(fEGGBinning)  = 11000;
    }
    delete arrEG;

    TString EGGBinning = GetEnvValue("EGGBinning",(TString)0);
    TObjArray *arrEGG = EGGBinning.Tokenize(" ");
    if(arrEGG->GetEntries() == 3) {
        get<0>(fEGGBinning) = ((TString)arrEGG->At(0)->GetName()).Atoi();
        get<1>(fEGGBinning) = ((TString)arrEGG->At(1)->GetName()).Atof();
        get<2>(fEGGBinning) = ((TString)arrEGG->At(2)->GetName()).Atof();
    }
    else {
        ERR_MESS << "Bad EGammaETot binning definition, needs 3 parameters, 8192 0 8192 used" << ENDL;
        get<0>(fEGGBinning)  = 11000;
        get<1>(fEGGBinning)  = 0;
        get<2>(fEGGBinning)  = 11000;
    }
    delete arrEGG;




    TString mode = GetEnvValue("Mode",(TString)0);

    fWorkingDir = GetEnvValue("WorkingDir",(TString)0);

    if(mode.Contains("TAG")) {
        TAG_Mult_ptr      = std::unique_ptr<TTreeReaderValue<Int_t>>(new TTreeReaderValue<Int_t>(fReader, "TAG_Mult"));
        TAG_QTot_ptr      = std::unique_ptr<TTreeReaderArray<Float_t>>(new TTreeReaderArray<Float_t>(fReader, "TAG_QTot"));
        TAG_QShort_ptr    = std::unique_ptr<TTreeReaderArray<Float_t>>(new TTreeReaderArray<Float_t>(fReader, "TAG_QShort"));
        TAG_IsFission_ptr = std::unique_ptr<TTreeReaderArray<Bool_t>>(new TTreeReaderArray<Bool_t>(fReader, "TAG_IsFission"));
        TAG_Id_ptr        = std::unique_ptr<TTreeReaderArray<Short_t>>(new TTreeReaderArray<Short_t>(fReader, "TAG_Id"));
        TAG_TS_ptr        = std::unique_ptr<TTreeReaderArray<Long64_t>>(new TTreeReaderArray<Long64_t>(fReader, "TAG_TS"));
    }
    else if(mode.Contains("LaBr3")){
        LaBr3_Mult_ptr      = std::unique_ptr<TTreeReaderValue<Int_t>>(new TTreeReaderValue<Int_t>(fReader, "LaBr3_Mult"));
        LaBr3_E_ptr    = std::unique_ptr<TTreeReaderArray<Float_t>>(new TTreeReaderArray<Float_t>(fReader, "LaBr3_E"));
        LaBr3_Id_ptr      = std::unique_ptr<TTreeReaderArray<Short_t>>(new TTreeReaderArray<Short_t>(fReader, "LaBr3_Id"));
        LaBr3_TS_ptr = std::unique_ptr<TTreeReaderArray<Long64_t>>(new TTreeReaderArray<Long64_t>(fReader, "LaBr3_TS"));
    }

    memset(fLast_AddMult,0,sizeof(fLast_AddMult));

    fPlotBasic = PlotOptions.Contains("Basic");
    fPlotMonitoring = PlotOptions.Contains("Monitoring");
    fPlotRates = PlotOptions.Contains("Rates");
    fPlotCalibrated = PlotOptions.Contains("Calibrated");
    fPlotRelAngles = PlotOptions.Contains("RelAngles");
    fPlotXTalk = PlotOptions.Contains("XTalk");
    fPlotDevHists = PlotOptions.Contains("Dev");
    fPlotOsloTSC = PlotOptions.Contains("OsloTSC");
    fUseSumThresholdForGG = PlotOptions.Contains("UseSumThresholdforGG");

    //fSumThresholdForGG = GetEnvValue("SumThresholdForGG",(TString)0);

    fEWindowEnergy = GetEnvValue("EWindowEnergy",(Int_t)0);
    fEWindowWidth = GetEnvValue("EWindowWidth",(Int_t)0);


    fReferenceRun = GetEnvValue("ReferenceRun",(Int_t)0);
    fFirstDate.Set(GetEnvValue("ReferenceTime",(TString)0));

    if(mode.Contains("TAG")) fUseTAG = true;
    if(mode.Contains("LaBr3")) {
        fUseLaBr3 = true;
    }

    fGGIsomers = GetEnvValue("GGIsomer",(TString)0);
    arr = fGGIsomers.Tokenize("&&");
    for(int i=0 ; i<arr->GetEntries() ; i++) {
        TString name = arr->At(i)->GetName();
        TObjArray *arr2 = name.Tokenize(" ");
        GGIsomerStruct astruct;
        if(arr2->GetEntries() == 2) {
            astruct.G1Min = ((TString)arr2->At(0)->GetName()).Atof();
            astruct.G1Max = ((TString)arr2->At(1)->GetName()).Atof();
            astruct.G2Min = ((TString)arr2->At(0)->GetName()).Atof();
            astruct.G2Max = ((TString)arr2->At(1)->GetName()).Atof();
        }
        else if(arr2->GetEntries() == 4) {
            astruct.G1Min = ((TString)arr2->At(0)->GetName()).Atof();
            astruct.G1Max = ((TString)arr2->At(1)->GetName()).Atof();
            astruct.G2Min = ((TString)arr2->At(2)->GetName()).Atof();
            astruct.G2Max = ((TString)arr2->At(3)->GetName()).Atof();
        }
        else {
            WARN_MESS << "bad GGIsomers definition in :" << name << " --> skipped" << ENDL;
            delete arr2;
            continue;
        }
        fListofGGIsomers.push_back(astruct);
        delete arr2;
    }
    delete arr;

    fGGEff_AngCorrCond = GetEnvValue("GGEff_AngCorrCond",(TString)0);
    arr = fGGEff_AngCorrCond.Tokenize("&&");
    for(int i=0 ; i<arr->GetEntries() ; i++) {
        TString name = arr->At(i)->GetName();
        TObjArray *arr2 = name.Tokenize(" ");
        if(arr2->GetEntries() !=4) {
            WARN_MESS << "bad GGEff_AngCorrCond definition in :" << name << " --> skipped" << ENDL;
            delete arr2;
            continue;
        }
        GGEffAngCorStruc astruct;
        astruct.EGate = ((TString)arr2->At(0)->GetName()).Atof();
        astruct.EWidth = ((TString)arr2->At(1)->GetName()).Atof();
        astruct.a2 = ((TString)arr2->At(2)->GetName()).Atof();
        astruct.a4 = ((TString)arr2->At(3)->GetName()).Atof();
        fListofGGEffAngCorrCond.push_back(astruct);
        delete arr2;
    }
    delete  arr;

    if(RelAnglesMode.Contains("clov",TString::ECaseCompare::kIgnoreCase)) fDoAngularUsingCloverAngles = true;

    if(fUseProof) {
        fArrayOfPositions = (TObjArray *)fInput->FindObject("PositionLUT");
    }

    if(fArrayOfPositions == nullptr) {
        fPlotRelAngles = false;
        ERR_MESS << "PositionLUT not found ==> Rel angles skipped" << ENDL;
    }
    else if(fDoAngularUsingCloverAngles) {
        TObjArray *temp = (TObjArray*)fArrayOfPositions->Clone();
        fArrayOfPositions->Clear();
        for(int i=0 ; i<temp->GetEntries()/4 ; i++) {
            TVector3 *vecclov = new TVector3;
            for(int j=0 ; j<4 ; j++) {
                TVector3 *vec2 = (TVector3*) temp->At(i*4+j);
                *vecclov += *vec2;
            }
            *vecclov *= 1./4.;
            fArrayOfPositions->Add(vecclov);
        }
    }

    if(fPlotRelAngles){
        TObjArray *arr = RelAnglesBinning.Tokenize(" ");
        if(arr->GetEntries() == 3) {
            get<0>(fRelAnglesBinning.first) = ((TString)arr->At(0)->GetName()).Atoi();
            get<1>(fRelAnglesBinning.first) = ((TString)arr->At(1)->GetName()).Atof();
            get<2>(fRelAnglesBinning.first) = ((TString)arr->At(2)->GetName()).Atof();

            get<0>(fRelAnglesBinning.second) = ((TString)arr->At(0)->GetName()).Atoi();
            get<1>(fRelAnglesBinning.second) = ((TString)arr->At(1)->GetName()).Atof();
            get<2>(fRelAnglesBinning.second) = ((TString)arr->At(2)->GetName()).Atof();
        }
        else if(arr->GetEntries() == 6) {
            get<0>(fRelAnglesBinning.first) = ((TString)arr->At(0)->GetName()).Atoi();
            get<1>(fRelAnglesBinning.first) = ((TString)arr->At(1)->GetName()).Atof();
            get<2>(fRelAnglesBinning.first) = ((TString)arr->At(2)->GetName()).Atof();
            get<0>(fRelAnglesBinning.second) = ((TString)arr->At(3)->GetName()).Atoi();
            get<1>(fRelAnglesBinning.second) = ((TString)arr->At(4)->GetName()).Atof();
            get<2>(fRelAnglesBinning.second) = ((TString)arr->At(5)->GetName()).Atof();
        }
        else {
            ERR_MESS << "Bad binning definition, needs 3 (symetric) or 6 (asymetric) parameters, 4096 0 4094 used" << ENDL;
            get<0>(fRelAnglesBinning.first)  = 4096;
            get<1>(fRelAnglesBinning.first)  = 0;
            get<2>(fRelAnglesBinning.first)  = 4096;
            get<0>(fRelAnglesBinning.second) = 4096;
            get<1>(fRelAnglesBinning.second) = 0;
            get<2>(fRelAnglesBinning.second) = 4096;
        }
        delete arr;

        // check for options
        if(RelAnglesNorm.Contains("normonly",TString::ECaseCompare::kIgnoreCase)) {
            fDoNormalizarion = true;
            fDoAngularCorrelations = false;
        }
        else if(RelAnglesNorm.Contains("norm",TString::ECaseCompare::kIgnoreCase)) {
            fDoNormalizarion = true;
            fDoAngularCorrelations = true;
        }
        else {
            fDoNormalizarion = false;
            fDoAngularCorrelations = true;
        }

        arr = RelAnglesCond.Tokenize(" ");

        for(int i=0 ; i<arr->GetEntries() ; i++){
            Float_t E,W;
            TString name = arr->At(i)->GetName();
            name.ToLower();
            if(name=="nocond"){
                E=0.;
                W=0.;
            }
            else if( name.Contains("_") ){
                TObjArray *arr2 = name.Tokenize("_");
                if(arr2->GetEntries()!=2){
                    delete arr2;
                    continue;
                }

                E = ((TString)arr2->First()->GetName()).Atof();
                W = ((TString)arr2->Last()->GetName()).Atof();
                delete arr2;
            }

            pair <Float_t, Float_t> apair(E,W);
            fRelAnglesGates.push_back(apair);

            map < Int_t, array < TH2*, 3 > > amap;
            fRelativeAnglesGG.push_back(amap);
            map < Int_t, array < TH2*, 3 > > amap_norm;
            fRelativeAnglesGG_Norm.push_back(amap_norm);

            pair <Int_t, Int_t> apair2(0,0);
            fRelAnglesNFired.push_back(apair2);
        }
        delete arr;

        if(fRelAnglesGates.size() == 0 ){
            fPlotRelAngles = false;
        }
    }

    /// ****************************** ///
    /// *** Define Basic Histograms *** ///
    /// ****************************** ///

    TString Folder;

    if(fPlotBasic){
        Folder = "Basic";

        hGlob_Mult       = BuildHist1D("Glob_Mult",       fMaxIds,0,fMaxIds,"Global multiplicity",       Folder   ,"TH1D");
        hMeanMultByEvt   = BuildHist1D("MeanMultByEvt",   12,0,12,          "Detector Type",               Folder   ,"TH1D");
        const char *Labels[11] = {"","FIIPS","IFIN","","FIPPS Back","FIPPS Front","FIPPS Side","IFIN AC","","TAG",""};
        for(int i=0 ; i<11 ; i++) hMeanMultByEvt->GetXaxis()->SetBinLabel(i+1,Labels[i]);

        Folder = "Basic/Ids";
        hADC_Values      = BuildHist1D("ADC_Values",       fMaxADC,0,fMaxADC,"ADC value",    Folder);
        hIds             = BuildHist1D("Ids",              fMaxIds,0,fMaxIds,"Channel Id",   Folder);
        hADC_vs_Ids      = BuildHist2D("ADC_vs_Ids",       fMaxIds,0,fMaxIds,"Channel Id",   fMaxADC,0,fMaxADC,"ADC value",  Folder);
        hIds_vs_Ids      = BuildHist2D("Ids_vs_Ids",       fMaxIds,0,fMaxIds,"Channel Id",   fMaxIds,0,fMaxIds,"Channel Id",  Folder);

        Folder = "Basic/RawTiming";

        hTiming_DT_to_Id0     = BuildHist2D("DT_to_Id0",   fTimeWindow/2,-fTimeWindow,fTimeWindow ,"Time between a detecor and Id==0 (ns)",                      fMaxIds,0,fMaxIds,"Detector Id", Folder);
        hTiming_DT_vsEnergyTAG  = BuildHist2D("hTiming_DT_vsEnergyTAG",   fTimeWindow/2,-fTimeWindow,fTimeWindow ,"Time between TAG and ref vs energy", 32000,0,32000,"Detector Id", Folder);
        hTiming_Window_DT     = BuildHist2D("Window_DT",   (fDTClover+fTimeWindow)/4,-fDTClover, fTimeWindow,"Clover Time related to the first event of the window (ns)",  fMaxNClov,0,fMaxNClov  ,"Clover Id",   Folder);
        if(fUseTAG) {
            hTiming_Clov_TAG  = BuildHist2D("DT_Clov-TAG", fTimeWindow/2,-fTimeWindow,fTimeWindow ,"Time between a Clover and the associated fission tag (ns)",  fMaxNClov,0,fMaxNClov  ,"Clover Id",   Folder);
            if(fNTAG>1) hTAG_vs_TAG       = BuildHist2D("TAG_vs_TAG",  1000,0,25000 ,"TAG 0",  1000,0,25000 ,"TAG 1",Folder);
        }
    }

    /// ************************************ ///
    /// *** Define Monitoring Histograms *** ///
    /// ************************************ ///

    Int_t   NumbOfDays = 40;
    Float_t BinWidthMonitorInSec = 600;
    Float_t BinWidthRateInSec = 60;

    fRateScaleFactor = (1./BinWidthRateInSec);

    Int_t NBinsMonitor = NumbOfDays*24*(3600/BinWidthMonitorInSec);
    Float_t BinMinMonitor=0;
    Float_t BinMaxMonitor=NumbOfDays*24*3600;

    Int_t NMaxRuns = 10000;

    if(fPlotMonitoring){
        Folder = "Monitoring";
        hEGammavsTime  = BuildHist2D("EGamma_vs_Time",  NBinsMonitor, BinMinMonitor, BinMaxMonitor ,"Date",        4000,0,8000,  "Energy (keV)",  Folder , true);

        if(fEWindowEnergy>0) {
            Folder = "Monitoring/EvsRun";
            for(int i=0 ; i<fMaxNClov*4 ; i++){
                hEGammaCalvsRun[i]    = BuildHist2D(Form("ECal-Run_%d",i),     NMaxRuns,fReferenceRun,fReferenceRun+ NMaxRuns      ,"Run number",  4*fEWindowWidth,fEWindowEnergy-fEWindowWidth,fEWindowEnergy+fEWindowWidth,  "ECal (keV)",Folder);
            }
        }

        if(fUseTAG){
            for(int i=0 ; i<fNTAG ; i++){
                hQTotvsTime[i]    = BuildHist2D(Form("QTot_vs_Time_%d",i)   ,NBinsMonitor, BinMinMonitor, BinMaxMonitor     ,"Date"       ,1000,0,30000, Form("QTot %d (channels)",i),Folder , true);
                hQTotvsRun[i]     = BuildHist2D(Form("QTot_vs_Run_%d",i)    ,NMaxRuns,fReferenceRun,fReferenceRun+ NMaxRuns ,"Run number" ,1000,0,30000, Form("QTot %d (channels)",i),Folder);
                histoQTot[i]      = BuildHist1D(Form("QTot_%d",i),32768,0,32768, Form("QTot %d (channels)",i),Folder,"TH1D");
            }
        }
    }

    /// ******************************* ///
    /// *** Define Rates Histograms *** ///
    /// ******************************* ///

    Int_t   NBinsRates = NumbOfDays*24*(3600/BinWidthRateInSec);
    Float_t BinMinRates=0;
    Float_t BinMaxRates=NumbOfDays*24*3600;

    if(fPlotRates){
        Folder = "Rates";
        hFIPPSRate         = BuildHist1D("FIPPSRate",       NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);
        hIFINRate          = BuildHist1D("IFINRate",        NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);
        hACRate            = BuildHist1D("AC_Rate",         NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);
        hADC45Rate         = BuildHist1D("ADC45Rate",       NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);

        if(fUseTAG){
            hTAGRate       = BuildHist1D("TAGRate",         NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);
            hFissionRate   = BuildHist1D("FissionRate",     NBinsRates, BinMinRates, BinMaxRates, "Date", Folder, "TH1F", true);
        }
    }
    /// ******************************************** ///
    /// *** Define Calibrated spectra Histograms *** ///
    /// ******************************************** ///

    if(fPlotCalibrated){
        // Float_t GGMax=4096;
        Int_t   GGNBins =   get<0>(fEGGBinning);
        Float_t GGMin   =   get<1>(fEGGBinning);
        Float_t GGMax   =   get<2>(fEGGBinning);

        Int_t   GGNBins_s=4096;
        Float_t GGMin_s=0;
        Float_t GGMax_s=2048;

        Int_t   GNBins  = get<0>(fEGBinning);
        Float_t GMin    = get<1>(fEGBinning);
        Float_t GMax    = get<2>(fEGBinning);


        if(fUseLaBr3){
            Folder = "Calibrated/LaBr3";
            for(int i=0 ; i<fNLaBr3 ; i++){
                hLaBr3_ESpectra[i]         = BuildHist1D(Form("LaBr3_%.2d",i)   ,  GNBins,GMin,GMax,   "Energy (keV)", Folder,"TH1D");
            }
        }

        Folder = "Calibrated/GeE";
        for(int i=0 ; i<fMaxNClov*4 ; i++){
            hGe_ESpectra[i]         = BuildHist1D(Form("GeE_%.2d",i)   ,  GNBins,GMin,GMax,   "Energy (keV)",    "Calibrated/CalE/All","TH1D");
            hGe_ESpectraMult[0][i]  = BuildHist1D(Form("GeE_M1_%.2d",i),  GNBins,GMin,GMax,   "Energy (keV)",    "Calibrated/CalE/M1","TH1D");
            hGe_ESpectraMult[1][i]  = BuildHist1D(Form("GeE_M2_%.2d",i),  GNBins,GMin,GMax,   "Energy (keV)",    "Calibrated/CalE/M2","TH1D");
        }

        for(int i=0 ; i<fMaxNClov ; i++)
            hSumE[i] = BuildHist1D(Form("SumE_%.2d",i),    GNBins,GMin,GMax,   "Energy (keV)",    "Calibrated/CalE/Sum","TH1D");

        for(int i=0 ; i<fMaxNClov ; i++){
            if(fPlotXTalk) {
                for(int ii=0 ; ii<4 ; ii++) {
                    for(int jj=ii+1 ; jj<4 ; jj++) {
                        hGe_EGGM2[i][ii][jj] = BuildHist2D(Form("XTalk_Clov%d_%d_%d",i,ii,jj),     GGNBins_s,GGMin_s,GGMax_s,  "E#gamma_{1} (keV)", GGNBins_s,GGMin_s,GGMax_s,   "E#gamma_{2} (keV)",   Form("Calibrated/XTalk/Clov%d",i));
                    }
                }
            }
        }

        hGe_EvsId = BuildHist2D("GeEvsId", GNBins,GMin,GMax,        "Energy (keV)",      fNGeDet,0,fNGeDet,   "Channel Id",        Folder);
        hGe_E     = BuildHist1D("GeE",    GNBins,GMin,GMax,       "Energy (keV)",      Folder,"TH1D");
        hGe_E_FIPPS = BuildHist1D("GeE_FIPPS",     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
        hGe_E_IFIN  = BuildHist1D("GeE_IFIN",     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
        hGe_E_GG  = BuildHist2D("GeE_GG",  GGNBins,GGMin,GGMax, "E#gamma_{1} (keV)",  GGNBins,GGMin,GGMax, "E#gamma_{2} (keV)", Folder);
        hGe_E_GG_FIPPS_FIPPS = BuildHist2D("GeE_GG_FIPPS_FIPPS",  GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{1} (keV)",  GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{2} (keV)", Folder);
        hGe_E_GG_IFIN_IFIN   = BuildHist2D("GeE_GG_IFIN_IFIN",    GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{1} (keV)",   GGNBins_s,GGMin_s,GGMax_s,"E#gamma_{2} (keV)", Folder);
        hGe_E_GG_FIPPS_IFIN  = BuildHist2D("GeE_GG_FIPPS_IFIN",  GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{1} (keV)",  GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{2} (keV)", Folder);
        hGe_E_GG_FIPPS_IFIN_sym  = BuildHist2D("GeE_GG_FIPPS_IFIN_sym",   GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{1} (keV)",   GGNBins_s,GGMin_s,GGMax_s, "E#gamma_{2} (keV)", Folder);

        for(auto i: fListofGGEffAngCorrCond) {
            array<TH1*,3> a_array;
            TString Name = Form("GeE_EffAngCorr_%.0f_%.0f_FIPPS_FIPPS",i.EGate,i.EWidth);
            a_array.at(0) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            Name = Form("GeE_EffAngCorr_%.0f_%.0f_FIPPS_IFIN",i.EGate,i.EWidth);
            a_array.at(1) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            Name = Form("GeE_EffAngCorr_%.0f_%.0f_IFIN_IFIN",i.EGate,i.EWidth);
            a_array.at(2) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            fListofGGEffAngCorrCondHists.push_back(a_array);
        }

        Folder = "Calibrated/Addback";
        for(int i=0 ; i<fMaxNClov ; i++){
            hAddE[i] = BuildHist1D(Form("AddE_%.2d",i),    GNBins,GMin,GMax,   "Energy (keV)",    Folder,"TH1D");
            hERej[i] = BuildHist1D(Form("EAC_Rej_%.2d",i), GNBins,GMin,GMax,   "Energy (keV)",    Folder,"TH1D");
        }

        for(auto i: fListofGGEffAngCorrCond) {
            array<TH1*,3> a_array_add;
            TString Name = Form("AddE_EffAngCorr_%.0f_%.0f_FIPPS_FIPPS",i.EGate,i.EWidth);
            a_array_add.at(0) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            Name = Form("AddE_EffAngCorr_%.0f_%.0f_FIPPS_IFIN",i.EGate,i.EWidth);
            a_array_add.at(1) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            Name = Form("AddE_EffAngCorr_%.0f_%.0f_IFIN_IFIN",i.EGate,i.EWidth);
            a_array_add.at(2) = BuildHist1D(Name,     GNBins,GMin,GMax,        "Energy (keV)",      Folder,"TH1D");
            fListofGGAddEffAngCorrCondHists.push_back(a_array_add);
        }

        Folder = "Calibrated";
        hEGammavsId  = BuildHist2D("EGammavsId",      GNBins,GMin,GMax,    "Energy (keV)", fMaxNClov,0,fMaxNClov,"Clover Id",  Folder);
        hEGamma      = BuildHist1D("EGamma",        GNBins,GMin,GMax,     "Energy (keV)", Folder,"TH1D");
        hEGamma_FIPPS      = BuildHist1D("hEGamma_FIPPS",        GNBins,GMin,GMax,     "Energy (keV)", Folder,"TH1D");
        hEGamma_IFIN      = BuildHist1D("EGamma_IFIN",        GNBins,GMin,GMax,     "Energy (keV)", Folder,"TH1D");
        hEGG         = BuildHist2D("EGG",           GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)",  GGNBins,GGMin,GGMax,  "E#gamma_{2} (keV)",   Folder);

        if(fUseSumThresholdForGG){
            hEGG_ESumThreshold = BuildHist2D("EGG_ESumThreshold", GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)",  GGNBins,GGMin,GGMax,  "E#gamma_{2} (keV)",   Folder);
        }
           
        hEGamma_ETot = BuildHist2D("EGamma_ETot",       get<0>(fEGammaETotBinning.first),get<1>(fEGammaETotBinning.first),get<2>(fEGammaETotBinning.first), "Calorimetric Energy (keV)" ,  get<0>(fEGammaETotBinning.second),get<1>(fEGammaETotBinning.second),get<2>(fEGammaETotBinning.second),  "E#gamma (keV)",   Folder);

        if(fPlotOsloTSC){
           hEGamma_ETot_M2 = BuildHist2D("EGamma_ETot_M2", get<0>(fEGammaETotBinning.first),get<1>(fEGammaETotBinning.first),get<2>(fEGammaETotBinning.first), "Calorimetric Energy (mult 2) (keV)" ,  get<0>(fEGammaETotBinning.second),get<1>(fEGammaETotBinning.second),get<2>(fEGammaETotBinning.second),  "E#gamma (keV)",   Folder);
        }

        if(fUseTAG){
            hEGamma_FTag = BuildHist1D("EGamma_FTag",  8000,0,8000,   "Energy (keV)", Folder,"TH1D");
            hEGamma_BTag = BuildHist1D("EGamma_BTag",  8000,0,8000,   "Energy (keV)", Folder,"TH1D");
            hEGG_FTag    = BuildHist2D("EGG_FTag",     GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)", GGNBins,GGMin,GGMax,   "E#gamma_{2} (keV)",   Folder);
            hEGG_BTag    = BuildHist2D("EGG_BTag",     GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)", GGNBins,GGMin,GGMax,   "E#gamma_{2} (keV)",   Folder);
            hEGamma_Qtot = BuildHist2D("hEGamma_hQTOT",     GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)", 32000,0,32000,   "QTAG ",   Folder);
            hEGamma_FTag_QTot = BuildHist2D("hEGamma_FTag_QTot",     GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)", 32000,0,32000,   "QTAG ",   Folder);
            hDeltaT_EGamma = BuildHist2D("hDeltaT_EGamma",  fTimeWindow/2,-fTimeWindow,fTimeWindow,  "DeltaT", 8000,0,8000,   "E#gamma_{1} (keV) ",   Folder);
     
        }

        Folder = "Calibrated/Isomers";
        for(auto &i: fListofGGIsomers) {
            TString Name = Form("GGIso_[%.0f-%.0f]_[%.0f-%.0f]",i.G1Min,i.G1Max,i.G2Min,i.G2Max);
            i.hist = BuildHist2D(Name,     GGNBins,GGMin,GGMax,  "E#gamma_{1} (keV)", GGNBins,GGMin,GGMax,   "E#gamma_{2} (keV)",   Folder);
        }
    }

    if(fPlotDevHists) {
        Folder = "Dev/GG_Clov_vs_Clov";

        for(int i=0 ; i<fMaxNClov ; i++){
            for(int j=i ; j<fMaxNClov ; j++){
                hMapOfClovbyClovGG[i][j] = BuildHist2D(Form("Clov%d_Clov%d",i,j),  4096,0,4096, Form("E#gamma Clov %d (keV)",i), 4096,0,4096, Form("E#gamma Clov %d (keV)",j), Folder);
            }
        }
    }

    /// ******************************************** ///
    /// *** Define Relative Angles Histograms *** ///
    /// ******************************************** ///

    if(fPlotRelAngles){
        TString modes[3] = {"FIPPS-FIPPS","FIPPS-IFIN","IFIN-IFIN"};

        hRelativeAnglesDist = new TH1F("RelAngles/RelAngleDist","RelAngleDist",182,0,182);
        fOutput->Add(hRelativeAnglesDist);

        fDummyMatrix = new TH2F("Dummy","Dummy",get<0>(fRelAnglesBinning.first),get<1>(fRelAnglesBinning.first),get<2>(fRelAnglesBinning.first),get<0>(fRelAnglesBinning.second),get<1>(fRelAnglesBinning.second),get<2>(fRelAnglesBinning.second));

        // Theta, N FIPPS-FIPPS, N FIPPS-IFIN, N IFIN-IFIN
        map <int, tuple<int,int,int> > combinations;
        map <int, TString > combinations_det;

        for(int idet=0 ; idet<fArrayOfPositions->GetEntries() ; idet++) {
            TVector3 *vec1 = (TVector3*) fArrayOfPositions->At(idet);
            if(!vec1 || vec1->Mag()==0.) continue;

            for(int idet2=0 ; idet2<=idet ; idet2++){
                TVector3 *vec2 = (TVector3*) fArrayOfPositions->At(idet2);
                if(!vec2 || vec2->Mag()==0.) continue;

                Int_t Theta = TMath::Nint(vec1->Angle(*vec2)*TMath::RadToDeg());
                if(Theta <= 28) continue;

                Int_t test = (((idet/4)>7) + ((idet2/4)>7));
                if(fDoAngularUsingCloverAngles) test = (((idet)>7) + ((idet2)>7));
                if(test==0) get<0>(combinations[Theta])++;
                if(test==1) get<1>(combinations[Theta])++;
                if(test==2) get<2>(combinations[Theta])++;

                hRelativeAnglesDist->Fill(Theta);

                TString comb;
                if(fDoAngularUsingCloverAngles) comb = Form("%.2d/%.2d  ",idet/4,idet2/4);
                else comb = Form("%.2d/%.2d  ",idet,idet2);
                if(!combinations_det[Theta].Contains(comb)) combinations_det[Theta] += comb;

                for(uint i=0 ; i<fRelAnglesGates.size() ; i++){
                    Float_t E = fRelAnglesGates[i].first;
                    Float_t W = fRelAnglesGates[i].second;

                    TString Name = Form("E%.1f_W%.1f",E,W);
                    if(E==0.)
                        Name = "NoCond";

                    TString Type = "Ge";
                    if(fDoAngularUsingCloverAngles) Type = "Clov";

                    if(fDoAngularCorrelations && (fRelativeAnglesGG[i].count(Theta) == 0 || fRelativeAnglesGG[i][Theta][test] == nullptr)) {
                        Folder = "RelAngles/" + Name + "/" + modes[test] + "/";
                        fRelativeAnglesGG[i][Theta][test] = new TH2F(Form("%s/GxG_%s_Theta_%d_%s_%s",Folder.Data(),Type.Data(),Theta,modes[test].Data(),Name.Data()),Form("GxG_%s_Theta_%d__%s_%s",Type.Data(),Theta,modes[test].Data(),Name.Data()),get<0>(fRelAnglesBinning.first),get<1>(fRelAnglesBinning.first),get<2>(fRelAnglesBinning.first),get<0>(fRelAnglesBinning.second),get<1>(fRelAnglesBinning.second),get<2>(fRelAnglesBinning.second));
                        fOutput->Add(fRelativeAnglesGG[i][Theta][test]);
                    }
                    if(fDoNormalizarion && (fRelativeAnglesGG_Norm[i].count(Theta) == 0 || fRelativeAnglesGG_Norm[i][Theta][test] == nullptr)) {
                        Folder = "RelAngles/" + Name + "/" + modes[test] + "/";
                        fRelativeAnglesGG_Norm[i][Theta][test] = new TH2F(Form("%s/GxG_%s_Theta_%d_%s_%s_Norm",Folder.Data(),Type.Data(),Theta,modes[test].Data(),Name.Data()),Form("GxG_%s_Theta_%d__%s_%s_Norm",Type.Data(),Theta,modes[test].Data(),Name.Data()),get<0>(fRelAnglesBinning.first),get<1>(fRelAnglesBinning.first),get<2>(fRelAnglesBinning.first),get<0>(fRelAnglesBinning.second),get<1>(fRelAnglesBinning.second),get<2>(fRelAnglesBinning.second));
                        fOutput->Add(fRelativeAnglesGG_Norm[i][Theta][test]);
                    }
                }
            }
        }

        for(auto i : combinations) {
            cout<<"Theta: "<<i.first<<" FIPPS-FIPPS:"<<get<0>(i.second)<<" FIPPS-IFIN:"<<get<1>(i.second)<<" IFIN-IFIN:"<<get<2>(i.second)<<endl;
        }
        cout<<"Combinations per angle:"<<endl;
        for(auto i : combinations_det) {
            cout<<"Theta="<<setw(3)<<i.first<<"=> "<<i.second<<endl;
        }

        if(fDoNormalizarion) {
            fLastPos_BinX.resize(NLastVecs);
            fLastPos_BinY.resize(NLastVecs);

            for(auto &i: fLastPos_BinX) {
                i.resize(get<0>(fRelAnglesBinning.first));
            }
            for(auto &i: fLastPos_BinY) {
                i.resize(get<0>(fRelAnglesBinning.second));
            }
        }
    }
}

TH1 *FIPPSSelector::BuildHist1D(TString HistName, Int_t NBins, Float_t Min, Float_t Max, TString AxisName, TString FolderName, TString Type, Bool_t Timing)
{
    TString Name;
    TH1 *hist = 0x0;

    Name = FolderName + "/" + HistName;

    if(Type=="TH1F")
        hist = new TH1F(Name,Name,NBins,Min,Max);
    else if(Type=="TH1I")
        hist = new TH1I(Name,Name,NBins,Min,Max);
    else if(Type=="TH1D")
        hist = new TH1D(Name,Name,NBins,Min,Max);

    hist->GetXaxis()->SetTitle(AxisName);

    hist->SetTitle(HistName);

    if(Timing){
        hist->GetXaxis()->SetTimeOffset(fFirstDate.Convert());
        hist->GetXaxis()->SetTimeDisplay(1);

        if(HistName.Contains("Rate"))
            hist->GetYaxis()->SetTitle("Rate (Hz)");
    }

    fOutput->Add(hist);

    return hist;
}

TH2 *FIPPSSelector::BuildHist2D(TString HistName,
                                Int_t NBinsX, Float_t MinX, Float_t MaxX, TString AxisNameX,
                                Int_t NBinsY, Float_t MinY, Float_t MaxY, TString AxisNameY,
                                TString FolderName,
                                Bool_t Timing)
{
    TString Name;
    TH2F *hist = 0x0;
    Name = FolderName + "/" + HistName;

    hist = new TH2F(Name,Name,NBinsX,MinX,MaxX,NBinsY,MinY,MaxY);
    hist->GetXaxis()->SetTitle(AxisNameX);
    hist->GetYaxis()->SetTitle(AxisNameY);

    hist->SetTitle(HistName);

    if(Timing){
        hist->GetXaxis()->SetTimeOffset(fFirstDate.Convert());
        hist->GetXaxis()->SetTimeDisplay(1);
    }

    fOutput->Add(hist);

    return hist;
}


TH3 *FIPPSSelector::BuildHist3D(TString HistName,
                                Int_t NBinsX, Float_t MinX, Float_t MaxX, TString AxisNameX,
                                Int_t NBinsY, Float_t MinY, Float_t MaxY, TString AxisNameY,
                                Int_t NBinsZ, Float_t MinZ, Float_t MaxZ, TString AxisNameZ,
                                TString FolderName)
{
    TString Name;
    TH3F *hist = 0x0;
    Name = FolderName + "/" + HistName;
    hist = new TH3F(Name,Name,NBinsX,MinX,MaxX,NBinsY,MinY,MaxY,NBinsZ,MinZ,MaxZ);
    hist->GetXaxis()->SetTitle(AxisNameX);
    hist->GetYaxis()->SetTitle(AxisNameY);
    hist->GetZaxis()->SetTitle(AxisNameZ);

    fOutput->Add(hist);

    return hist;
}

void FIPPSSelector::SlaveTerminate()
{
}

void FIPPSSelector::Terminate()
{
    if(fUseProof == false){
        printf("\nSelector::Process() is done. Now writing output file.\n"); fflush(stdout);

        if(fPlotBasic)
            hMeanMultByEvt->Scale(1./hGlob_Mult->GetEntries());

        TObjArray *arr;
        for(int i=0 ; i<fOutput->GetEntries() ; i++){
            if(!(fOutput->At(i)->InheritsFrom("TH1") || fOutput->At(i)->InheritsFrom("TGraph")  || fOutput->At(i)->InheritsFrom("TF1")))
                continue;

            TString Name = fOutput->At(i)->GetName();
            arr = Name.Tokenize("/");
            TString FileName = arr->Last()->GetName();
            TString FolderName="";
            for(int j=0 ; j<arr->GetEntries()-1 ; j++)
                FolderName += Form("%s/",arr->At(j)->GetName());

            fFinalFile->mkdir(FolderName);
            fFinalFile->cd(FolderName);
            ((TNamed*)fOutput->At(i))->SetName(FileName);
            fOutput->At(i)->Write(fOutput->At(i)->GetName(),TObject::kSingleKey);
            delete arr;
        }

        if(fFinalFile->GetListOfKeys()->GetEntries()==0){
            gSystem->Unlink(fFinalFile->GetName());
            fFinalFile->Delete();
            fFinalFile = nullptr;
        }
        else{
            cout<<"\e[1;92mOutput FileName = "<<fFileNameOut<<"\e[0m"<<endl<<endl;
            fFinalFile->Close();
        }
    }

    if(fUseProof == false) delete progress;
}


Bool_t FIPPSSelector::Process(Long64_t entry)
{
    fReader.SetLocalEntry(entry);
    bool stop = false;
    // Check Run modification
    if(fReader.GetTree()->GetTree() != fCurrentTree) {
        fCurrentTree = fReader.GetTree()->GetTree();
        fCurrentRun = atoi(fCurrentTree->GetUserInfo()->FindObject("Run")->GetTitle());
    }

    if(fUseTAG) {
        TAG_Mult = TAG_Mult_ptr.get()->Get();
        TAG_QTot = (Float_t*)TAG_QTot_ptr.get();
        TAG_QShort = (Float_t*)TAG_QShort_ptr.get();
        TAG_IsFission = (Bool_t*)TAG_IsFission_ptr.get();
        TAG_Id = (Short_t*)TAG_Id_ptr.get();
        TAG_TS = (Long64_t*) TAG_TS_ptr.get();
    }
    if(fUseLaBr3) {
        LaBr3_Mult = LaBr3_Mult_ptr.get()->Get();
        LaBr3_E = (Float_t*)LaBr3_E_ptr.get();
        LaBr3_Id = (Short_t*)LaBr3_Id_ptr.get();
        LaBr3_TS = (Long64_t*) LaBr3_TS_ptr.get();
    }

    if(fPlotBasic) {
        Long64_t TS_Id0=0;

        for(int i=0 ; i<*multiplicity ; i++) {
            if(globid[i]==0) {
                TS_Id0 = timeStamp[i];
                break;
            }
        }

        // Analyse raw data
        hGlob_Mult->Fill(*multiplicity);

        for(int i=0 ; i<*multiplicity ; i++) {
            hADC_Values->Fill(adc[i]);
            hIds->Fill(globid[i]);
            hADC_vs_Ids->Fill(globid[i],adc[i]);
            for(int j=0 ; j<i ; j++){
                hIds_vs_Ids->Fill(globid[i],globid[j]);
                hIds_vs_Ids->Fill(globid[j],globid[i]);
            }
            if(TS_Id0 && !(timeStamp[i]==TS_Id0 && globid[i]==0)) hTiming_DT_to_Id0->Fill(timeStamp[i]-TS_Id0,globid[i]);
            
            
            if(globid[i]>=96) hTiming_DT_vsEnergyTAG->Fill(timeStamp[i]-TS_Id0, TAG_QTot[0]);

            if(adc[i]==45 && fPlotRates) hADC45Rate->Fill(*GlobTime,fRateScaleFactor);
        }
    }

    Int_t FissionId=-1;

    // Ana TAG
    if(fUseTAG && *TAG_Mult) {
        for(int itag=0 ; itag<*TAG_Mult ; itag++) {
            if(TAG_IsFission[itag]) {
                FissionId = itag;
                break;
            }
        }
        if(fPlotBasic) {
            hMeanMultByEvt->Fill(9,*TAG_Mult);
            if(*TAG_Mult==2 && TAG_Id[0] != TAG_Id[1] && hTAG_vs_TAG ) {
                hTAG_vs_TAG->Fill(TAG_QTot[0],TAG_QTot[1]);
                hTAG_vs_TAG->Fill(TAG_QTot[1],TAG_QTot[0]);
            }
        }
        if(fPlotMonitoring) {
            hQTotvsTime[(Int_t)TAG_Id[0]]->Fill(*GlobTime,TAG_QTot[0]);
            hQTotvsRun[(Int_t)TAG_Id[0]]->Fill(fCurrentRun,TAG_QTot[0]);
            histoQTot[(Int_t)TAG_Id[0]]->Fill(TAG_QTot[0]);
        }
        if(fPlotRates) {
            hTAGRate->Fill(*GlobTime,fRateScaleFactor);
            if(FissionId>=0) hFissionRate->Fill(*GlobTime,fRateScaleFactor);
        }
    }
    
    //AnaLaBr3
    if(fUseLaBr3 && *LaBr3_Mult) {
        for(int iLaBr=0 ; iLaBr<*LaBr3_Mult ; iLaBr++) {
            Int_t LaBr3Id = LaBr3_Id_ptr.get()->At(iLaBr);
            if(fPlotCalibrated){
                hLaBr3_ESpectra[LaBr3Id]->Fill(LaBr3_E_ptr.get()->At(iLaBr)); 
            }
        };
    }

    // Ana Ge
    for(int ige=0 ; ige<*Ge_Mult ; ige++) {
        Int_t GeId = Ge_ClovId[ige]*4+Ge_Id[ige];
        if(fPlotBasic) {
            if(Ge_Type[ige] == kFIPPS) hMeanMultByEvt->Fill(1);
            else if(Ge_Type[ige] == kIFIN) hMeanMultByEvt->Fill(2);
        }
        if(fPlotMonitoring) {
            if(fEWindowEnergy>0 && TMath::Abs(Ge_E[ige]-fEWindowEnergy)<fEWindowWidth) hEGammaCalvsRun[GeId]->Fill(fCurrentRun,Ge_E[ige]);
        }
        if(fPlotRates) {
            if(Ge_Type[ige] == kFIPPS) hFIPPSRate->Fill(*GlobTime,fRateScaleFactor);
            if(Ge_Type[ige] == kIFIN) hIFINRate->Fill(*GlobTime,fRateScaleFactor);
        }
        if(fPlotCalibrated){
            hGe_ESpectra[GeId]->Fill(Ge_E[ige]);
            if(Ge_ClovMult[ige]==1) hGe_ESpectraMult[0][GeId]->Fill(Ge_E[ige]);
            else if(Ge_ClovMult[ige]==2) hGe_ESpectraMult[1][GeId]->Fill(Ge_E[ige]);

            if(fPlotXTalk && Ge_ClovMult[ige]==2 && Ge_ClovId[ige] == Ge_ClovId[ige+1] && TMath::Abs(Ge_TS[ige]-Ge_TS[ige+1])<fDTClover) {
                hGe_EGGM2[Ge_ClovId[ige]][Ge_Id[ige]][Ge_Id[ige+1]]->Fill(Ge_E[ige],Ge_E[ige+1]);
                hGe_EGGM2[Ge_ClovId[ige]][Ge_Id[ige+1]][Ge_Id[ige]]->Fill(Ge_E[ige+1],Ge_E[ige]);
            }

            hGe_EvsId->Fill(Ge_E[ige],GeId);
            hGe_E->Fill(Ge_E[ige]);
            if(((Int_t)Ge_Type[ige]) == kFIPPS) hGe_E_FIPPS->Fill(Ge_E[ige]);
            if(((Int_t)Ge_Type[ige]) == kIFIN) hGe_E_IFIN->Fill(Ge_E[ige]);
            hSumE[Ge_ClovId[ige]]->Fill(Ge_E[ige]);

            for(int ige2=0 ; ige2<ige ; ige2++) {
                hGe_E_GG->Fill(Ge_E[ige],Ge_E[ige2]);
                hGe_E_GG->Fill(Ge_E[ige2],Ge_E[ige]);

                if(((Int_t)Ge_Type[ige])==kFIPPS && ((Int_t)Ge_Type[ige2])==kFIPPS) {
                    hGe_E_GG_FIPPS_FIPPS->Fill(Ge_E[ige],Ge_E[ige2]);
                    hGe_E_GG_FIPPS_FIPPS->Fill(Ge_E[ige2],Ge_E[ige]);
                }
                if(((Int_t)Ge_Type[ige])==kIFIN && ((Int_t)Ge_Type[ige2])==kIFIN) {
                    hGe_E_GG_IFIN_IFIN->Fill(Ge_E[ige],Ge_E[ige2]);
                    hGe_E_GG_IFIN_IFIN->Fill(Ge_E[ige2],Ge_E[ige]);
                }
                if(((Int_t)Ge_Type[ige])==kFIPPS && ((Int_t)Ge_Type[ige2])==kIFIN) {
                    hGe_E_GG_FIPPS_IFIN->Fill(Ge_E[ige],Ge_E[ige2]);

                    hGe_E_GG_FIPPS_IFIN_sym->Fill(Ge_E[ige],Ge_E[ige2]);
                    hGe_E_GG_FIPPS_IFIN_sym->Fill(Ge_E[ige2],Ge_E[ige]);
                }
                if(((Int_t)Ge_Type[ige])==kIFIN && ((Int_t)Ge_Type[ige2])==kFIPPS) {
                    hGe_E_GG_FIPPS_IFIN->Fill(Ge_E[ige2],Ge_E[ige]);

                    hGe_E_GG_FIPPS_IFIN_sym->Fill(Ge_E[ige],Ge_E[ige2]);
                    hGe_E_GG_FIPPS_IFIN_sym->Fill(Ge_E[ige2],Ge_E[ige]);
                }
            }
            for(size_t i=0 ; i<fListofGGEffAngCorrCond.size() ; i++) {
                GGEffAngCorStruc *astruc = &fListofGGEffAngCorrCond.at(i);
                if( abs(Ge_E[ige]-astruc->EGate)<astruc->EWidth) {
                    for(int ige2=0 ; ige2<*Ge_Mult ; ige2++) {
                        if(ige==ige2) continue;
                        TVector3 vec1(Ge_X[ige],Ge_Y[ige],Ge_Z[ige]);
                        TVector3 vec2(Ge_X[ige2],Ge_Y[ige2],Ge_Z[ige2]);
                        if(vec2.Mag()==0. || vec1.Mag()==0.) continue;
                        Double_t CosTheta = TMath::Cos(vec1.Angle(vec2));
                        Double_t Weight = 1.+astruc->a2*0.5*(3.*CosTheta*CosTheta-1)+astruc->a4*(1./8.)*(35*pow(CosTheta,4)-30.*CosTheta*CosTheta+3.);

                        if(((Int_t)Ge_Type[ige])==kFIPPS && ((Int_t)Ge_Type[ige2])==kFIPPS) {
                            fListofGGEffAngCorrCondHists.at(i).at(0)->Fill(Ge_E[ige2],1./Weight);
                        }
                        if(((Int_t)Ge_Type[ige])==kIFIN && ((Int_t)Ge_Type[ige2])==kIFIN) {
                            fListofGGEffAngCorrCondHists.at(i).at(2)->Fill(Ge_E[ige2],1./Weight);
                        }
                        if(((Int_t)Ge_Type[ige])==kFIPPS && ((Int_t)Ge_Type[ige2])==kIFIN) {
                            fListofGGEffAngCorrCondHists.at(i).at(1)->Fill(Ge_E[ige2],1./Weight);
                        }
                        if(((Int_t)Ge_Type[ige])==kIFIN && ((Int_t)Ge_Type[ige2])==kFIPPS) {
                            fListofGGEffAngCorrCondHists.at(i).at(1)->Fill(Ge_E[ige2],1./Weight);
                        }
                    }
                }
            }
        }
    };

    // Ana AC
    for(int iac=0 ; iac<*AC_Mult ; iac++) {
        if(fPlotBasic) {
            if(AC_Type[iac] == kFIPPS) hMeanMultByEvt->Fill(4+AC_Id[iac]);
            else if(AC_Type[iac] == kIFIN) hMeanMultByEvt->Fill(7);
        }
        if(fPlotRates) {
            hACRate->Fill(*GlobTime,fRateScaleFactor);
        }
    }

    // Ana Clover
    if(fPlotRelAngles) {
        for(uint i=0 ; i<fRelAnglesNFired.size() ; i++){
            fRelAnglesNFired[i].first = 0;
            fRelAnglesNFired[i].second = 0;
        }
        for(int iclov=0 ; iclov<*Add_Mult ; iclov++) {
            for(uint ig=0 ; ig<fRelAnglesGates.size() ; ig++){
                Float_t E = fRelAnglesGates[ig].first;
                Float_t W = fRelAnglesGates[ig].second;

                if( E==0. || (TMath::Abs(Add_E[iclov]-E) < W )){
                    fRelAnglesNFired[ig].first++;
                    fRelAnglesNFired[ig].second=Add_ClovId[iclov];
                }
            }
        }
    }

    // calc calorimetric energy
    Float_t ECalor = 0;
    Float_t ECalorM2 = 0;
    for(int iclov=0 ; iclov<*Add_Mult ; iclov++) {
        ECalor += Add_E[iclov];
	
	    if(*Add_Mult == 2) 
            ECalorM2 += Add_E[iclov];

        if(fPlotCalibrated) {
            //efficiency plots
            for(size_t i=0 ; i<fListofGGEffAngCorrCond.size() ; i++) {
                GGEffAngCorStruc *astruc = &fListofGGEffAngCorrCond.at(i);
                if( abs(Add_E[iclov]-astruc->EGate)<astruc->EWidth) {
                    for(int iclov2=0 ; iclov2<*Add_Mult ; iclov2++) {
                        if(iclov==iclov2) continue;
                        TVector3 vec1(Add_X[iclov],Add_Y[iclov],Add_Z[iclov]);
                        TVector3 vec2(Add_X[iclov2],Add_Y[iclov2],Add_Z[iclov2]);
                        if(vec2.Mag()==0. || vec1.Mag()==0.) continue;
                        Double_t CosTheta = TMath::Cos(vec1.Angle(vec2));
                        Double_t Weight = 1.+astruc->a2*0.5*(3.*CosTheta*CosTheta-1)+astruc->a4*(1./8.)*(35*pow(CosTheta,4)-30.*CosTheta*CosTheta+3.);

                        Double_t Mode=1.;
                        int id1 = Add_ClovId[iclov];
                        int id2 = Add_ClovId[iclov2];
                        if(id1<8 && id2<8) Mode=0.;
                        else if(id1>=8 && id2>=8) Mode=2.;

                        fListofGGAddEffAngCorrCondHists.at(i).at(Mode)->Fill(Add_E[iclov2],1./Weight);
                    }
                }
            }
        }
    }

    for(int iclov=0 ; iclov<*Add_Mult ; iclov++) {
        if(fPlotBasic) {
            if(Add_TS[iclov] != timeStamp[0]) hTiming_Window_DT->Fill(Add_TS[iclov]-timeStamp[0],Add_ClovId[iclov]);
            if(fUseTAG && *TAG_Mult) hTiming_Clov_TAG->Fill(Add_TS[iclov]-TAG_TS[0],Add_ClovId[iclov]);

            if(fUseTAG && *TAG_Mult && Add_E[iclov]>5.) {
            Long64_t TS_Id0=0, TS_Id96=0;

            for(int i=0 ; i<*multiplicity ; i++) {
                if(globid[i]==0) {
                    TS_Id0 = timeStamp[i];
                }
                if(globid[i]>=96){
                    TS_Id96 = timeStamp[i];
                }
            }
            hDeltaT_EGamma->Fill(TS_Id96-TS_Id0,Add_E[iclov]);
            }
            
        }
        if(fPlotMonitoring) {
            if(Add_E[iclov]>0.) hEGammavsTime->Fill(*GlobTime, Add_E[iclov]);
        }
        if(fPlotDevHists) {
            if(Add_E[iclov]>0.) {
                for(int iclov2=0 ; iclov2<iclov ; iclov2++) {
                    if(Add_E[iclov2]>0.) {
                        Int_t id1, id2;
                        if(Add_ClovId[iclov]<Add_ClovId[iclov2]) {
                            id1 = iclov;
                            id2 = iclov2;
                        }
                        else {
                            id1 = iclov2;
                            id2 = iclov;
                        }

                        Int_t Clov1 = Add_ClovId[id1];
                        Int_t Clov2 = Add_ClovId[id2];

                        Int_t E1 = Add_E[id1];
                        Int_t E2 = Add_E[id2];

                        hMapOfClovbyClovGG[Clov1][Clov2]->Fill(E1,E2);
                    }
                }
            }
        }
        if(fPlotCalibrated || fPlotRelAngles) {
            if(Add_E[iclov]>0.) {

                if(fPlotCalibrated) {
                    hAddE[Add_ClovId[iclov]]->Fill(Add_E[iclov]);

                    if(Add_ClovId[iclov] < 8) 
			        hEGamma_FIPPS->Fill(Add_E[iclov]);
                    if(Add_ClovId[iclov] >= 8) 
                        hEGamma_IFIN->Fill(Add_E[iclov]);

                    hEGammavsId->Fill(Add_E[iclov],Add_ClovId[iclov]);
                    hEGamma->Fill(Add_E[iclov]);
                    hEGamma_ETot->Fill(ECalor,Add_E[iclov]);
                    
                    if(fUseTAG) {
                        if(FissionId>=0) hEGamma_FTag->Fill(Add_E[iclov]);
                        else hEGamma_BTag->Fill(Add_E[iclov]);

                        hEGamma_Qtot -> Fill (Add_E[iclov], TAG_QTot[0]);
                    if(FissionId>=0)   hEGamma_FTag_QTot -> Fill(Add_E[iclov], TAG_QTot[0]);

                    //hDeltaT_EGamma->Fill(timeStamp[i]-TS_Id0,Add_E[iclov]);
                    }
                }
                for(int iclov2=0 ; iclov2<iclov ; iclov2++) {
                    if(Add_E[iclov2]>0.) {
                        if(fPlotCalibrated) {
                            hEGG->Fill(Add_E[iclov],Add_E[iclov2]);
                            hEGG->Fill(Add_E[iclov2],Add_E[iclov]);


                            if(fUseSumThresholdForGG && (Add_E[iclov]+Add_E[iclov2]) > 3000){
                                hEGG_ESumThreshold->Fill(Add_E[iclov],Add_E[iclov2]);
                                hEGG_ESumThreshold->Fill(Add_E[iclov2],Add_E[iclov]);
                            }

                            for(auto i: fListofGGIsomers) {
                                if(i.IsValid(Add_TS[iclov]-Add_TS[0],Add_TS[iclov2]-Add_TS[0])) {
                                    i.hist->Fill(Add_E[iclov],Add_E[iclov2]);
                                    i.hist->Fill(Add_E[iclov2],Add_E[iclov]);
                                }
                            }

                            if(fUseTAG && FissionId>=0) {
                                hEGG_FTag->Fill(Add_E[iclov],Add_E[iclov2]);
                                hEGG_FTag->Fill(Add_E[iclov2],Add_E[iclov]);
                            }
                            else if(fUseTAG && FissionId<0) {
                                hEGG_BTag->Fill(Add_E[iclov],Add_E[iclov2]);
                                hEGG_BTag->Fill(Add_E[iclov2],Add_E[iclov]);
                            }

                            if(fPlotOsloTSC){
                                if (*Add_Mult == 2) {
                                    hEGamma_ETot_M2->Fill(ECalorM2,Add_E[iclov]);
                                    hEGamma_ETot_M2->Fill(ECalorM2,Add_E[iclov2]);
                                }
                            }
                        }

                        //Relative Angle Calculations
                        if(fDoAngularCorrelations) {

                            // Couple d'énergie E1 E2
                            Float_t E1 = Add_E[iclov];
                            Float_t E2 = Add_E[iclov2];

                            // Bins correspondant au couple E1 sur X, E2 sur Y
                            Int_t BinE1X = fDummyMatrix->GetXaxis()->FindBin(E1)-1;
                            Int_t BinE2Y = fDummyMatrix->GetYaxis()->FindBin(E2)-1;

                            // Bins correspondant au couple E2 sur X, E1 sur Y
                            Int_t BinE2X = fDummyMatrix->GetXaxis()->FindBin(E2)-1;
                            Int_t BinE1Y = fDummyMatrix->GetYaxis()->FindBin(E1)-1;

                            bool ok = false;
                            bool cond1 = false;
                            bool cond2 = false;

                            // Check if we are within the window
                            if(( (BinE1X >=0) && (BinE1X < fDummyMatrix->GetNbinsX())) && ((BinE2Y >=0) && (BinE2Y < fDummyMatrix->GetNbinsY()) )) cond1 = true;
                            if(( (BinE1Y >=0) && (BinE1Y < fDummyMatrix->GetNbinsY())) && ((BinE2X >=0) && (BinE2X < fDummyMatrix->GetNbinsX()) )) cond2 = true;

                            ok = (cond1 || cond2);

                            if(ok) {
                                TVector3 vec1;
                                TVector3 vec2;

                                vec1.SetXYZ(Add_X[iclov],Add_Y[iclov],Add_Z[iclov]);
                                vec2.SetXYZ(Add_X[iclov2],Add_Y[iclov2],Add_Z[iclov2]);

                                // recuperation des positions vec1 (E1), et vec2(E2)
                                if(fDoAngularUsingCloverAngles && max(Add_ClovId[iclov],Add_ClovId[iclov2]) < fArrayOfPositions->GetEntries()) {
                                    vec1 = *((TVector3*)fArrayOfPositions->At(Add_ClovId[iclov]));
                                    vec2 = *((TVector3*)fArrayOfPositions->At(Add_ClovId[iclov2]));
                                }

                                // on supprime les positions de detecteur non definis
                                if(vec1.Mag()==0 || vec2.Mag()==0) continue;

                                if((fUseTAG && FissionId>=0) || !fUseTAG) {
                                    //boucle sur les gates
                                    for(uint i=0 ; i<fRelAnglesGates.size() ; i++) {
                                        bool doangcor = false;
                                        // si gate valide plus d'une fois, on rempli tous les couples
                                        if(fRelAnglesNFired[i].first>1) doangcor = true;
                                        // si gate valide une seule fois, on rempli que les couples qui ne contiennent pas l'energie de la gate
                                        if(fRelAnglesNFired[i].first==1 && (Add_ClovId[iclov] != fRelAnglesNFired[i].second && Add_ClovId[iclov2] != fRelAnglesNFired[i].second) ) doangcor = true;

                                        if(doangcor) {
                                            Int_t Theta = TMath::Nint(vec1.Angle(vec2)*TMath::RadToDeg());
                                            if(fRelativeAnglesGG[i].count(Theta)) {
                                                Double_t Mode=1.;
                                                int id1 = Add_ClovId[iclov];
                                                int id2 = Add_ClovId[iclov2];
                                                if(id1<8 && id2<8) Mode=0.;
                                                else if(id1>=8 && id2>=8) Mode=2.;

                                                if(cond1) fRelativeAnglesGG[i][Theta][Mode]->Fill(E1,E2);
                                                if(cond2) fRelativeAnglesGG[i][Theta][Mode]->Fill(E2,E1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                if(fPlotCalibrated) hERej[Add_ClovId[iclov]]->Fill(Add_ERej[iclov]);
            }
        }
    }
    if(fPlotRelAngles && fDoNormalizarion && NFilled<NLastVecs) {
        fLast_AddMult[NFilled] = *Add_Mult;
        for(int iclov=0 ; iclov<fLast_AddMult[NFilled] ; iclov++) {
            fLast_Add_E[NFilled][iclov] = Add_E[iclov];
            fLast_Add_X[NFilled][iclov] = Add_X[iclov];
            fLast_Add_Y[NFilled][iclov] = Add_Y[iclov];
            fLast_Add_Z[NFilled][iclov] = Add_Z[iclov];
            fLast_Add_TS[NFilled][iclov] = Add_TS[iclov];
            fLast_Add_ClovId[NFilled][iclov] = Add_ClovId[iclov];
        }
        NFilled++;
    }
    if(fPlotRelAngles && fDoNormalizarion && NFilled==NLastVecs) {
        for(int ipast1=0 ; ipast1<NLastVecs ; ipast1++) {
            for(int ipast2=0 ; ipast2<ipast1 ; ipast2++) {

                for(int iclov1=0 ; iclov1<fLast_AddMult[ipast1] ; iclov1++) {
                    Float_t E1 = fLast_Add_E[ipast1][iclov1];
                    if(E1==0.) continue;
                    for(int iclov2=0 ; iclov2<fLast_AddMult[ipast2] ; iclov2++) {
                        Float_t E2 = fLast_Add_E[ipast2][iclov2];
                        if(E2==0.) continue;

                        // Bins correspondant au couple E1 sur X, E2 sur Y
                        Int_t BinE1X = fDummyMatrix->GetXaxis()->FindBin(E1)-1;
                        Int_t BinE2Y = fDummyMatrix->GetYaxis()->FindBin(E2)-1;

                        // Bins correspondant au couple E2 sur X, E1 sur Y
                        Int_t BinE2X = fDummyMatrix->GetXaxis()->FindBin(E2)-1;
                        Int_t BinE1Y = fDummyMatrix->GetYaxis()->FindBin(E1)-1;

                        bool ok = false;
                        bool cond1 = false;
                        bool cond2 = false;

                        // Check if we are within the window
                        if(( (BinE1X >=0) && (BinE1X < fDummyMatrix->GetNbinsX())) && ((BinE2Y >=0) && (BinE2Y < fDummyMatrix->GetNbinsY()) )) cond1 = true;
                        if(( (BinE1Y >=0) && (BinE1Y < fDummyMatrix->GetNbinsY())) && ((BinE2X >=0) && (BinE2X < fDummyMatrix->GetNbinsX()) )) cond2 = true;

                        ok = (cond1 || cond2);

                        if(ok) {
                            TVector3 vec1;
                            TVector3 vec2;

                            vec1.SetXYZ(fLast_Add_X[ipast1][iclov1],fLast_Add_Y[ipast1][iclov1],fLast_Add_Z[ipast1][iclov1]);
                            vec2.SetXYZ(fLast_Add_X[ipast2][iclov2],fLast_Add_Y[ipast2][iclov2],fLast_Add_Z[ipast2][iclov2]);

                            // recuperation des positions vec1 (E1), et vec2(E2)
                            if(fDoAngularUsingCloverAngles && max(fLast_Add_ClovId[ipast1][iclov1],fLast_Add_ClovId[ipast2][iclov2]) < fArrayOfPositions->GetEntries()) {
                                vec1 = *((TVector3*)fArrayOfPositions->At(fLast_Add_ClovId[ipast1][iclov1]));
                                vec2 = *((TVector3*)fArrayOfPositions->At(fLast_Add_ClovId[ipast2][iclov2]));
                            }

                            // on supprime les positions de detecteur non definis
                            if(vec1.Mag()==0 || vec2.Mag()==0) continue;

                            if((fUseTAG && FissionId>=0) || !fUseTAG) {
                                //boucle sur les gates
                                for(uint i=0 ; i<fRelAnglesGates.size() ; i++) {
                                    bool doangcor = false;
                                    // si gate valide plus d'une fois, on rempli tous les couples
                                    if(fRelAnglesNFired[i].first>1) doangcor = true;
                                    // si gate valide une seule fois, on rempli que les couples qui ne contiennent pas l'energie de la gate
                                    if(fRelAnglesNFired[i].first==1 && (fLast_Add_ClovId[ipast1][iclov1] != fRelAnglesNFired[i].second && fLast_Add_ClovId[ipast2][iclov2] != fRelAnglesNFired[i].second) ) doangcor = true;

                                    if(doangcor) {
                                        Int_t Theta = TMath::Nint(vec1.Angle(vec2)*TMath::RadToDeg());
                                        if(fRelativeAnglesGG_Norm[i].count(Theta)) {
                                            Double_t Mode=1.;
                                            int id1 = fLast_Add_ClovId[ipast1][iclov1];
                                            int id2 = fLast_Add_ClovId[ipast2][iclov2];
                                            if(id1<8 && id2<8) Mode=0.;
                                            else if(id1>=8 && id2>=8) Mode=2.;

                                            if(cond1) fRelativeAnglesGG_Norm[i][Theta][Mode]->Fill(E1,E2);
                                            if(cond2) fRelativeAnglesGG_Norm[i][Theta][Mode]->Fill(E2,E1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        NFilled=0;
    }

    if(stop) {
        cout<<endl;
        cin.get();
    }

    if(fUseProof == false) ++*progress;

    return kTRUE;
}
