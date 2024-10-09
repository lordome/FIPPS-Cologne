#include "LSTClovBuilder.h"

#include "Riostream.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TSystemDirectory.h"

#include "progressbar.h"
#include "CrateBoard.h"
#include "Event.h"
#include "EventArray.h"

using namespace lstdpp128;

int main(int argc, char** argv)
{
    if(argc<=1) {
        PrintHelp();
    }

    if(getenv("LSTClovBuilderConfFile") != nullptr && gSystem->IsFileInIncludePath(getenv("LSTClovBuilderConfFile"))) {
        fConfFileName = getenv("LSTClovBuilderConfFile");
    }

    bool doconfFile = false;
    bool doFileName = false;
    bool doNLoops = false;

    for(int i=1 ; i<argc ; i++) {
        TString par = argv[i];
        par.ToLower();
        if(par.Contains("conf",TString::kIgnoreCase) && !doconfFile) {doconfFile = true; continue;}
        else if(par.Contains("ow",TString::kIgnoreCase)) {fOverWrite = true; continue;}
        else if(par.Contains("nc",TString::kIgnoreCase)) {fNoCompress = true; continue;}
        else if(par.Contains("of",TString::kIgnoreCase)) {doFileName = true; continue;}
        else if(par.Contains("print",TString::kIgnoreCase)) {fPrintEvts = true; continue;}
        else if(par.Contains("batch",TString::kIgnoreCase)) {fBatchMode = true; continue;}
        else if(par.Contains("loops",TString::kIgnoreCase)) {doNLoops = true; continue;}
        else if(par.Contains("help",TString::kIgnoreCase)) {PrintHelp();}
        if(doFileName){
            TString FileName = argv[i];
            fOutputFileNameInit = FileName;
            doFileName = false;
            continue;
        }
        if(doNLoops){
            fNLoops = (Int_t)((TString)argv[i]).Atof();
            doNLoops=false;
            continue;
        }
        if(doconfFile){
            TString FileName = argv[i];
            if(gSystem->IsFileInIncludePath(FileName)) fConfFileName = FileName;
            else {
                ERR_MESS << "Configuration file: " << FileName << " not found" << endl;
                return 1;
            }
            doconfFile = false;
        }
    }

    // to define the extrenal variale
    ListModeContext listModeContext;

    if(ReadConf()) exit(EXIT_FAILURE);
    if(ReadLUT()) exit(EXIT_FAILURE);

    bool dorunlist = false;

    for(int i=1 ; i<argc ; i++) {
        TString par = argv[i];
        par.ToLower();
        if(par.Contains("run",TString::kIgnoreCase)) {dorunlist = true; continue;}
        if(dorunlist){
            TString runlist = argv[i];
            TObjArray *arr = runlist.Tokenize("-");
            fFirstRun = ((TString)arr->First()->GetName()).Atoi();
            if(runlist.EndsWith("-*")) {
                fInfiniteRuns = true;
                fLastRun = fFirstRun;
            }
            else if(runlist.Contains("-*")) {
                Int_t NRuns = ((TString)arr->Last()->GetName()).ReplaceAll("*","").Atoi();
                fLastRun = fFirstRun + NRuns-1;
            }
            else {
                fLastRun = ((TString)arr->Last()->GetName()).Atoi();
            }
            for(int i=fFirstRun ; i<=fLastRun ; i++) {

                TString FileName =  Form("%s/%.06d.root",fRawDataDir.Data(),fCurrentRun);               
                
                cout << FileName << endl;
                if(gSystem->IsFileInIncludePath(FileName))
                    fRuns.push_back(i);
            }
            dorunlist = false;
        }
    }

    cout<<endl;
    std::cout<<"\e[0;3;92m";
    cout<<"****************************"<<endl;
    cout<<"*  Analyse runs to treat   *"<<endl;
    cout<<"****************************"<<endl;
    cout<<ENDL;

    if(fOverWrite) WARN_MESS << " overwrite option used, existing files will be recreate !!!" << ENDL;
    if(fPrintEvts) WARN_MESS << " Print Evts flag activated !!!" << ENDL;

    if(fOutputFileNameInit!="" && fRuns.size()>1) {
        ERR_MESS << "User defined name cannot be used for many runs ==> EXIT" << ENDL;
        exit(EXIT_FAILURE);
    }

    if(fRuns.size()) {
        if(fInfiniteRuns) {
            INFO_MESS<<" Every runs from run " << fRuns[0] << " will be treated"<<ENDL;
        }
        else if(fRuns.size() == 1)
            INFO_MESS<<fRuns.size()<<" run to be treated: "<<fRuns[0]<<ENDL;
        else
            INFO_MESS<<fRuns.size()<<" runs to be treated, from "<<fRuns[0]<<" to "<<fRuns[fRuns.size()-1]<<ENDL;
    }
    else {
        WARN_MESS<<"No runs to treat... do you really need my services ?"<<ENDL;
        return 1;
    }
    cout<<endl;

    if(fBatchMode == false) {
        INFO_MESS << "Press a key to continue" << ENDL;
        cin.get();
    }

    gROOT->SetBatch();
    if(fInfiniteRuns) {
        fLastRun = GetLastRun();
        fCurrentRun = fFirstRun;
        while(true) {
            if(fCurrentRun<=fLastRun) {
                TString FileName = Form("%s/%.06d.root",fRawDataDir.Data(),fCurrentRun);   

                /*int test = */system(Form("ls %s > /dev/null",FileName.Data()));
                if(gSystem->IsFileInIncludePath(FileName)) {
                    INFO_MESS<<Form("Run: %.06d (remaining %d)",fCurrentRun,fLastRun-fCurrentRun)<<ENDL;
                    ReadFile();
                }
            }
            else {
                Int_t test = GetLastRun();
                if(test>fLastRun) {
                    fLastRun = test;
                    continue;
                }
                else break;
            }
            fCurrentRun++;
        }
    }
    else {
        for(uint i=0 ; i<fRuns.size() ; i++) {
            fCurrentRun = fRuns[i];
            INFO_MESS<<Form("Run: %.06d (remaining %d)",fCurrentRun,(Int_t)(fRuns.size()-i))<<ENDL;
            ReadFile();
        }
    }

    return 0;
}

void PrintHelp()
{
    cout<<endl;
    std::cout<<"\e[0;3;92m";
    cout << "*******************************"<<endl;
    cout << "** How to use LSTClovBuilder **"<<endl;
    cout << "********************************"<<endl;

    cout << "./LSTClovBuilder -options -runs" << endl;
    cout << "Available options: " << endl;
    cout << " --> -conf ConfFileName  ==> Define the Configuration file (if not specified, takes the LSTClovBuilder.conf in the current directory)" << endl;
    cout << " --> -ow                 ==> Overwrite the existing output file" << endl << endl;
    cout << " --> -of Name            ==> Define specific output file name" << endl << endl;
    cout << " --> -print              ==> Dump the events in the terminal" << endl << endl;
    cout << " --> -batch              ==> Batch mode, no need to press a key to validate the process" << endl << endl;
    cout << " --> -nc                 ==> No output file compression" << endl << endl;
    cout << " --> -nloops             ==> number of loops (4096 evts) to process" << endl << endl;

    cout << "Different ways to give the range of runs to treat:" << endl;
    cout << "-run 1234                ==> Will treat only the run nb 1234" << endl;
    cout << "-run 1234-2345           ==> Will treat all runs between 1234 and 2345" << endl;
    cout << "-run 1234-*              ==> Infinite loop from run 1234" << endl;
    cout << "-run 1234-*50            ==> Will treat 50 runs, starting with run 1234" << endl;
    cout<<ENDL;
    exit(0);
}

Int_t GetLastRun()
{
    TSystemDirectory dir(fRawDataDir,fRawDataDir);
    TList *l = dir.GetListOfFiles();
    while(l==nullptr) {
        WARN_MESS<<"Disk is busy, retry in few seconds"<<ENDL;
        sleep(2);
        l = dir.GetListOfFiles();
    }
    l->Sort();
    Int_t entry = l->GetEntries()-1;
    while(entry>=0) {
        TString Name = l->At(entry)->GetName();
        if(Name.EndsWith(".root")) {
            TObjArray *arr = Name.Tokenize("/");
            Name = arr->Last()->GetName();
            delete arr;
            Name.ReplaceAll("*.root","");
            return Name.Atoi();
        }
        entry--;
    }
    return -1;
}

Int_t ReadConf()
{
    cout<<endl;
    std::cout<<"\e[0;3;92m";
    cout<<"************************************"<<endl;
    cout<<"* Reading LSTClovBuilder conf File *"<<endl;
    cout<<"************************************"<<endl;
    cout<<ENDL;

    TITLE_MESS << "Conf File:" <<"\e[0;3;94m " << fConfFileName << ENDL << endl;

    std::ifstream FileConf;
    FileConf.open(fConfFileName.Data(),std::ifstream::in);

    if(!FileConf) {
        ERR_MESS<<"No Configuration file ==> Exit"<<ENDL;
        return 1;
    }

    string line;
    TString Buffer;

    while(FileConf) {
        getline(FileConf,line);
        Buffer = line;
        Buffer.ReplaceAll("\t"," ");
        if(Buffer.BeginsWith("#")) continue;
        TObjArray *arr=Buffer.Tokenize(" ");

        if(Buffer.BeginsWith("LUT_FILE",TString::kIgnoreCase)) {
            fLUTFile = arr->At(1)->GetName();
        }
        if(Buffer.BeginsWith("RawDataDir",TString::kIgnoreCase)) {
            fRawDataDir = arr->At(1)->GetName();
        }
        if(Buffer.BeginsWith("OutputDir",TString::kIgnoreCase)) {
            fOutputDir = arr->At(1)->GetName();
            gSystem->mkdir(fOutputDir,true);
        }
        if(Buffer.BeginsWith("DTClover",TString::kIgnoreCase)) {
            fDTClover = ((TString)arr->At(1)->GetName()).Atoi();
        }
        if(Buffer.BeginsWith("TimeBinning",TString::kIgnoreCase)) {
            fTimeBinning = ((TString)arr->At(1)->GetName()).Atoi();
        }
        if(Buffer.BeginsWith("DTAC",TString::kIgnoreCase)) {
            fDTACMin = ((TString)arr->At(1)->GetName()).Atoi();
            fDTACMax = ((TString)arr->At(2)->GetName()).Atoi();
        }
         if(Buffer.BeginsWith("idref",TString::kIgnoreCase)) {
            fIdRef = ((TString)arr->At(1)->GetName()).Atoi();
        }
        if(Buffer.BeginsWith("TimeAlignmentWindow",TString::kIgnoreCase)) {
            fTimeAlignmentWindow = ((TString)arr->At(1)->GetName()).Atoi();
        }
        if(Buffer.BeginsWith("PlotM",TString::kIgnoreCase)) {
            MultPlots = Buffer.Copy().ReplaceAll("PlotM ","");
            if(MultPlots.Contains("1")) fPlotM1 = true;
            if(MultPlots.Contains("2")) fPlotM2 = true;
            if(MultPlots.Contains("3")) fPlotM3 = true;
            if(MultPlots.Contains("4")) fPlotM4 = true;
        }
        delete arr;
    }

    if(fDTACMin == 0 && fDTACMax ==0) {
        fDTACMin = -fDTClover;
        fDTACMax = fDTClover;
    }

    TITLE_MESS <<"PATHS:" << ENDL;
    INFO_MESS<<"Raw data folder: "<<fRawDataDir<<ENDL;
    INFO_MESS<<"Output folder  : "<<fOutputDir<<ENDL;
    INFO_MESS<<"LUT File       : "<<fLUTFile<<ENDL;

    TITLE_MESS <<"Parameters:" << ENDL;

    INFO_MESS<<"NLoops  : "<<fNLoops<<ENDL;
    INFO_MESS<<"DTClover: "<<fDTClover<<ENDL;
    INFO_MESS<<"DTAC (min,max): "<<fDTACMin<<" "<<fDTACMax<<ENDL;
    INFO_MESS<<"Time window for time alignment: "<<fTimeAlignmentWindow<<ENDL;
    INFO_MESS<<"Clover multiplicity plots: " << MultPlots << ENDL;

    return 0;
}

Int_t ReadLUT()
{
    cout<<endl;
    std::cout<<"\e[0;3;92m";
    cout<<"****************************"<<endl;
    cout<<"*     Reading LUT File     *"<<endl;
    cout<<"****************************"<<endl;
    cout<<ENDL;

    TString Buffer;
    string line;

    std::ifstream FileConf;
    FileConf.open(fLUTFile, std::ifstream::in);

    if(!FileConf){
        ERR_MESS<<fLUTFile<<" not found ==> Exit"<<ENDL;
        return 1;
    }

    Int_t linenb=0;
    while(FileConf) {
        getline(FileConf,line);
        linenb++;
        Buffer = line;

        if(Buffer.BeginsWith("#") || Buffer.Length()==0 ) continue;

        TObjArray *loa=Buffer.ReplaceAll("\t"," ").Tokenize(" ");
        if(loa->GetEntries()==0) continue;
        if(loa->GetEntries() != 8) {
            ERR_MESS<<"Error in LUT.txt at line "<<linenb<<" ; skipped"<<ENDL;
            delete loa;
            continue;
        }

        Int_t adcval = ((TString)loa->At(0)->GetName()).Atoi();
        Int_t type = ((TString)loa->At(1)->GetName()).Atoi();
        Int_t clover = ((TString)loa->At(2)->GetName()).Atoi();
        Int_t detid = ((TString)loa->At(3)->GetName()).Atoi();
        Int_t globid = ((TString)loa->At(4)->GetName()).Atoi();
        Int_t TimeOffset  = ((TString)loa->At(5)->GetName()).Atoi();
        Int_t RangeMin  = ((TString)loa->At(6)->GetName()).Atoi();
        Int_t RangeMax  = ((TString)loa->At(7)->GetName()).Atoi();
        delete loa;

        DetDef *det = new DetDef;
        if(fDetectors.count(adcval)!=0 && det->Type>0) {
            WARN_MESS << "ADC "<<adcval<<" already seen in LUT !!! "<<ENDL;
            continue;
        }

        fDetectors[adcval] = det;

        det->Type = type;
        det->adc = adcval;
        det->CloverId = clover;
        det->Id = detid;
        det->GlobId = globid;
        det->TimeOffset = TimeOffset;
        det->RangeMin = RangeMin;
        det->RangeMax = RangeMax;

        if(det->IsFIPPS()) fNFIPPS++;
        else if(det->IsIFIN()) fNIFIN++;
        else if(det->IsTAG()) {
            fNTAG++;
            fUseTAG = true;
            det->Id = clover;
        }
        else if(det->IsLaBr3()) {
            fNLaBr3++;
            fUseLaBr3 = true;
            det->Id = clover;
        }
        else if(det->IsTAC()) {
            fNTAC++;
            fUseTAC = true;
            det->Id = clover;
        }
        else if(det->IsSiPM()) {
            fNSiPM++;
            fUseSiPM = true;
            det->Id = clover;
        }
        if(globid>fNIds) fNIds = globid;
    }

    map < Int_t, DetDef* >::iterator itr;
    for (itr = fDetectors.begin(); itr != fDetectors.end(); ++itr) {
        DetDef *det = itr->second;

        if(fDetectorsFromIds.count(det->GlobId)!=0 && det->Type>0) {
            WARN_MESS << "Id "<<det->GlobId<<" already seen in LUT !!! "<<ENDL;
            continue;
        }
        if(det->Type<=0) continue;

        fDetectorsFromIds[det->GlobId] = det;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].adc = det->adc;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].DetId = det->Id;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].GlobId = det->GlobId;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].DetType = det->Type;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].rangemax = det->RangeMax;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].rangemin = det->RangeMin;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].RefId = det->CloverId;
        fGlobHeader.LookUpTable.LUTEntry[det->GlobId].TimeOffset = det->TimeOffset;

        if(det->GlobId == 0) {
            fTimeOffsetId0 = det->TimeOffset;
            if(fTimeOffsetId0>0) WARN_MESS << "Det Glob id 0 (time reference) has a non null time offset of " << fTimeOffsetId0 << "Please check that this is normal" << ENDL;
        }
    }

    INFO_MESS<<fDetectors.size()<<" entries Read in the LUT"<<ENDL;

    FileConf.close();

    return 0;
}

TString GetRunDate(TString date)
{
    TObjArray *loa=date.ReplaceAll("\t"," ").Tokenize(" ");
    if(loa->GetEntries() != 2){
        cout<<"Error in RunsDates file ==> skipped"<<endl;
        delete loa;
        return "Undefined";
    }

    TString thedate = loa->At(0)->GetName();
    TObjArray *loa2=thedate.Tokenize("-");
    Int_t year = 2000 + ((TString)loa2->At(2)->GetName()).Atoi();
    TString months = loa2->At(1)->GetName();
    Int_t month=0;
    if(months == "Jan") month = 1;
    if(months == "Feb") month = 2;
    if(months == "Mar") month = 3;
    if(months == "Apr") month = 4;
    if(months == "May") month = 5;
    if(months == "Jun") month = 6;
    if(months == "Jul") month = 7;
    if(months == "Aug") month = 8;
    if(months == "Sep") month = 9;
    if(months == "Oct") month = 10;
    if(months == "Nov") month = 11;
    if(months == "Dec") month = 12;

    Int_t day = ((TString)loa2->At(0)->GetName()).Atoi();
    delete loa2;

    TString thetime = loa->At(1)->GetName();
    loa2=thetime.Tokenize(":");
    Int_t hour = ((TString)loa2->At(0)->GetName()).Atoi();
    Int_t min = ((TString)loa2->At(1)->GetName()).Atoi();
    Int_t sec = ((TString)loa2->At(2)->GetName()).Atoi();
    delete loa2;
    delete loa;

    TDatime datime(year,month,day,hour,min,sec);
    return datime.AsSQLString();
}


void ReadFile()
{
    gSystem->mkdir(Form("%s",fOutputDir.Data()),true);
    if(fOutputFileNameInit == "") fOutputFileName = Form("%s/%.06d.sort",fOutputDir.Data(),fCurrentRun);
    else fOutputFileName = fOutputFileNameInit;

    TString OutputFileNameCompressed = fOutputFileName + ".lz4";

    if(fOverWrite==false && (gSystem->IsFileInIncludePath(fOutputFileName) || gSystem->IsFileInIncludePath(OutputFileNameCompressed))) {
        if(gSystem->IsFileInIncludePath(OutputFileNameCompressed)) fOutputFileName = OutputFileNameCompressed;
        INFO_MESS << Form("Output file %s already existing ==> run %.06d skipped",fOutputFileName.Data(),fCurrentRun) << ENDL;
        return;
    }

    fOutputFile.open(fOutputFileName.Data());

    // Extract run date
    TString MetaFileName = Form("%s%.06d",fRawDataDir.Data(),fCurrentRun);
    ifstream file(MetaFileName);
    TString Buffer; string line;
    getline(file,line);getline(file,line);getline(file,line);
    Buffer=line;
    Buffer.Remove(0,14); // remove the user name
    Buffer.ReplaceAll("\t","");
    file.close();

    fCurrentDate = GetRunDate(Buffer);

    INFO_MESS << "Run date: " << fCurrentDate << ENDL;

    strcpy(fGlobHeader.date,fCurrentDate.Data());
    fGlobHeader.Run = fCurrentRun;
    fGlobHeader.DTClover = fDTClover;
    fGlobHeader.DTACMin = fDTACMin;
    fGlobHeader.DTACMax = fDTACMax;

    // Write header
    fOutputFile.write((char*)&fGlobHeader, sizeof(GlobHeader));
    //    fGlobHeader.Print();

    // Read input file
    TString FileName = Form("%s/%.06d.root",fRawDataDir.Data(),fCurrentRun); 

    cout <<  FileName << endl;
    fCurrentFile.open(FileName.Data());//, ios::in | ios::binary);
    if(!fCurrentFile){
        WARN_MESS<<FileName<<"not found"<<ENDL;
        return;
    }

    InitHistos();

    TFile *f = new TFile(FileName);

    TTree* t=(TTree*)f->Get("Data_R");

    UShort_t channel, board, energy;
    ULong64_t timestamp;
    UInt_t flags;
    t->SetBranchAddress("Channel",&channel);
    t->SetBranchAddress("Timestamp",&timestamp);
    t->SetBranchAddress("Board",&board);
    t->SetBranchAddress("Flags",&flags);
    t->SetBranchAddress("Energy",&energy);

    int32_t fsize =  1 << 12;
    int32_t fsizeInBytes = fsize * sizeof(int32_t);

    // Read the 4 fixed elements of the header.
    uint32_t nbEvents = t->GetEntries();

    // nbEvents = 100;
    uint32_t nbBoards = 2;
    uint32_t totalNumberOfBoards;
    int32_t *boardHeader;


    // Read the board structure.
    boardHeader = new int32_t[nbBoards];
    // fCurrentFile.read(reinterpret_cast<char *>(boardHeader), nbBoards * sizeof(uint32_t));

    // Read the list mode context.
    // totalNumberOfBoards = readListModeContextFromBinary(listModeContext, version, timeBase, nbEvents, nbBoards, (uint32_t*)boardHeader);
    totalNumberOfBoards = 2;

    INFO_MESS << listModeContext << ENDL;

    // creating the events that will be ordered
    Int_t BSize = nbEvents+1;
    //    Int_t BSize = nbEvents/10.;
    EventArray events(BSize);

    int n=0;

    progressbar progress(nbEvents);
    int fprocessedevents[11]; // 0: all, 1: FIPPS Ge, 2: IFIN Ge, 3 IFIN AC, 4: TAG, 5: FIPPS AC Back, 6: FIPPS AC Front, 7:  FIPPS AC Side, 8:  LaBr3, 9:  TAC, 10:  SiPM
    memset(fprocessedevents,0,sizeof(fprocessedevents));

    int fprocessedbadevents=0;
    int fprocessedThrRejections[11];
    int fprocessedPileUp[11];
    int fprocessednotinsertedevents=0;

    memset(fNGoodEvents,0,sizeof(fNGoodEvents));
    memset(fprocessedThrRejections,0,sizeof(fprocessedThrRejections));
    memset(fprocessedPileUp,0,sizeof(fprocessedPileUp));

       
        for (int index = 0; index<nbEvents; index++){          

            //if(index > 200 ) break;
            t->GetEntry(index);
            
            // reading event
            Event event;
           // bool testread = readEvent(event, fbuffer + index);

            event.crate = 0;
            event.channel = channel;
            event.board = board;
            event.adc = 1*event.board+event.channel;

            event.data[0] = energy;
        
            event.setTime(timestamp/1000);
            //flag for the moment is ignored. 

            int64_t eventts = event.time();

            int32_t adc = 8*event.board+event.channel;
            DetDef *det = fDetectors[adc];

            //cout << "ADC: " << adc << " channel: " << channel << " board " << board << endl;

            event.adc = adc;
            event.setTime(eventts+det->TimeOffset);//+det->TimeOffset);
            event.is_in_trigger = true ;
            EventArray::InsertResult result = events.insertEvent(event);

            if (result.process) {
                // do the process
                EventArray::array_type::iterator center = BuildEvents(result.begin, result.end, false);
                // update the processed center
                events.setProcessedCenter(center);
            }

        
            fprocessedevents[0]++;
            ++progress;
        }

        // if(fNLoops>0 && n==fNLoops) break;
        // n++;
    // }

    events.pack();

    // in case of init state, only process the beginning of the array
    // if (std::distance(events.center(), events.array().end()) < 0) {
    //     BuildEvents(events.array().begin(), events.array().end(), true);
    // } else {
        // processing in two steps, because the processor cannot manage blocks
        // with size greater than half the size of the events array
        EventArray::array_type::iterator center = BuildEvents(events.array().begin(), events.center(), false);

        // update the processed center
        events.setProcessedCenter(center);

        BuildEvents(center, events.array().end(), true);
    // }

    cout<<endl;

    INFO_MESS << Form("%5.5g",(double)fprocessedevents[0]) << " processed events " << ENDL;
    if(fUseTAC)  cout << "  --> TAC     : " << Form("%5.5g",(double)fprocessedevents[9]) << endl;
    if(fUseSiPM) cout << "  --> SiPM    : " << Form("%5.5g",(double)fprocessedevents[10]) << endl;

    if(fprocessedThrRejections[0]) {
        cout << "  --> Thresholds rejections: " << setw(8) << Form("%5.5g",(double)fprocessedThrRejections[0]);
        cout << "(";
        if(fprocessedThrRejections[6]) cout << Form(" TAC: %8.3g", ((double)fprocessedThrRejections[6]));
        if(fprocessedThrRejections[7]) cout << Form(" SiPM: %8.3g", ((double)fprocessedThrRejections[7]));
        cout << ")" << endl;
    }
    if(fprocessedPileUp[0]) {
        cout << "  --> Pile-up rejections:    " << setw(8) << Form("%5.5g",(double)fprocessedPileUp[0]);
        cout << "(";
        if(fprocessedPileUp[6]) cout << Form(" TAC: %8.3g", ((double)fprocessedPileUp[6]));
        if(fprocessedPileUp[7]) cout << Form(" SiPM: %8.3g", ((double)fprocessedPileUp[7]));
        cout << ")" << endl;
    }
    if(fprocessedThrRejections[1]) WARN_MESS << "  --> Unused ADC present in the dataflow: " << Form("%5.5g",(double)fprocessedThrRejections[1]) << ENDL;

    cout<<endl;
    INFO_MESS << Form("%5.5g",(double)fNGoodEvents[0]) << " clover events: " << ENDL;

    if(fUseTAC)   cout << "  --> TAC     : " << Form("%5.5g",(double)fNGoodEvents[9]) << endl;
    if(fUseSiPM)  cout << "  --> SiPM    : " << Form("%5.5g",(double)fNGoodEvents[10]) << endl;

    if(fprocessedbadevents) WARN_MESS << setw(10) << fprocessedevents << " processed bad events " << ENDL;
    if(fprocessednotinsertedevents) WARN_MESS << setw(10) << fprocessednotinsertedevents << " processed not insertable events " << ENDL;

    fCurrentFile.close();
    fOutputFile.close();

    WriteHistos();

    if(fNoCompress) {
        WARN_MESS << "lz4 compressor not used, output file will not be compressed" << ENDL;
    }
    else if(system("which lz4 >/dev/null")) {
        WARN_MESS << "lz4 compressor not found, output file will not be compressed" << ENDL;
    }
    else {
        INFO_MESS << "Compressing the output file..." << ENDL;
        int test = system(Form("lz4 -f %s %s.lz4",fOutputFileName.Data(),fOutputFileName.Data()));
        test = system(Form("rm -f %s",fOutputFileName.Data()));
        fOutputFileName = Form("%s.lz4",fOutputFileName.Data());
    }

    INFO_MESS << "Output Sorted file written: " << fOutputFileName << ENDL ;
}


template<typename Iterator>
Iterator BuildEvents(Iterator first, Iterator last, bool lastBlock)
{
    cout << endl << " ==> Building sorted events in terms of Clover events:"<<endl;
    progressbar progress(distance(first,last));

    Iterator position = first;
    Event *CurrentEvt = &(*(position));

    Iterator LastInWindowPosition = first;
    Event *LastInWindow = &(*(LastInWindowPosition));

    DetDef *det = nullptr;

    while (distance(position,last)>0) {
        CurrentEvt = &(*(position));
        det = fDetectors.at(CurrentEvt->adc);

        Long64_t CurrentTS = (Long64_t)CurrentEvt->time();

        if(fPrintEvts) {
            cout << left << "ADC: " << setw(4) << det->adc << " Glob Id: " << setw(4) <<  det->GlobId << " Type: " << setw(13) << fTypeString[det->Type] << " TS: " << setw(8) << CurrentTS << " energy: " << setw(6) << CurrentEvt->data[0] << " To be treated: " << CurrentEvt->is_in_trigger << endl;
        }

        if(det->IsClover()) {
            fGeHists[det->GlobId]->Fill(CurrentEvt->data[0]);
        }
        else if(det->IsAC()) fACHists[det->GlobId]->Fill(CurrentEvt->data[0]);
        else if(fUseTAG && det->IsTAG()) {
            fTAG_QTot[det->Id]->Fill(CurrentEvt->data[0]);
            fTAG_QShort[det->Id]->Fill(CurrentEvt->data[1]);
        }
        else if(fUseLaBr3 && det->IsLaBr3()) {
            fLaBr3_E[det->Id]->Fill(CurrentEvt->data[0]);
        }
        else if(fUseTAC && det->IsTAC()) {
            fTAC_T[det->Id]->Fill(CurrentEvt->data[0]);
        }
        else if(fUseSiPM && det->IsSiPM()) {
            fSiPM_T[det->Id]->Fill(CurrentEvt->fine_timestamp);
            fSiPM_QTot[det->Id]->Fill(CurrentEvt->data[0]);
            fSiPM_QShort[det->Id]->Fill(CurrentEvt->data[1]);
        }

        if(fLastTS>0)
            fDeltaTS->Fill(CurrentTS-fLastTS);
        if(fLastTS_per_id[det->GlobId]>0)
            fDeltaTSvsId->Fill(CurrentTS-fLastTS_per_id[det->GlobId],det->GlobId);
        fLastTS = CurrentTS;
        fLastTS_per_id[det->GlobId] = CurrentTS;
        if(det->IsClover()){
            if(fLastGeTS>0) fDeltaTS_Ge->Fill(CurrentTS-fLastGeTS);
            fLastGeTS = CurrentTS;
        }
        else if(det->IsAC()) {
            if(fLastACTS>0) fDeltaTS_AC->Fill(CurrentTS-fLastACTS);
            fLastACTS = CurrentTS;
        }
        else if(det->IsTAG()) {
            if(fLastTAGTS>0) fDeltaTS_TAG->Fill(CurrentTS-fLastTAGTS);
            fLastTAGTS = CurrentTS;
        }
        else if(det->IsLaBr3()) {
            if(fLastLaBr3TS>0) fDeltaTS_LaBr3->Fill(CurrentTS-fLastLaBr3TS);
            fLastLaBr3TS = CurrentTS;
        }
        else if(det->IsTAC()) {
            if(fLastTACTS>0) fDeltaTS_TAC->Fill(CurrentTS-fLastTACTS);
            fLastTACTS = CurrentTS;
        }
        else if(det->IsSiPM()) {
            if(fLastSiPMTS>0) fDeltaTS_TSiPM->Fill(CurrentTS-fLastSiPMTS);
            fLastSiPMTS = CurrentTS;
        }

        if(det->GlobId !=  fIdRef){
            if(CurrentTS-fLastTSId0 < fTimeAlignmentWindow){
                fTSDiffvsIdAlign->Fill(CurrentTS-fLastTSId0,det->GlobId);
                fTSDiffvsIdRaw->Fill((CurrentTS-det->TimeOffset)-(fLastTSId0-fTimeOffsetId0),det->GlobId);
            }
            else {
                pair<Int_t, ULong64_t > a_pair(det->GlobId,CurrentTS);
                fSavedTS.push_back(a_pair);

                while(CurrentTS-fSavedTS.front().second>fTimeAlignmentWindow) fSavedTS.erase(fSavedTS.begin());
            }
        }
        else {
            while(fSavedTS.size()) {
                if(CurrentTS-fSavedTS.front().second<fTimeAlignmentWindow){
                    fTSDiffvsIdAlign->Fill(-1.*(CurrentTS-fSavedTS.front().second),fSavedTS.front().first);
                    fTSDiffvsIdRaw->Fill(-1.*((CurrentTS-det->TimeOffset)-(fSavedTS.front().second-fDetectorsFromIds[fSavedTS.front().first]->TimeOffset)),fSavedTS.front().first);
                }
                fSavedTS.erase(fSavedTS.begin());
            };
            fLastTSId0 = CurrentTS;
        }

        if(CurrentEvt->is_in_trigger == false) {
            ++progress;
            position++;
            continue;
        }

        while(CurrentTS - LastInWindow->time() > fDTClover) {
            LastInWindowPosition++;
            LastInWindow = &(*(LastInWindowPosition));
        }

        if(det->IsClover()) {
            DetHeader header;
            Clov aclov;
            if(det->IsFIPPS()) header.type = kFIPPS;
            else header.type = kIFIN;
            aclov.ClovID = det->CloverId;

            Iterator testpos = LastInWindowPosition;
            Event *testevt = &(*(testpos));
            DetDef *testdet = fDetectors.at(testevt->adc);

            TString GeIds="";
            TString ACIds="";
            Int_t ACMultInClov = 0;
            vector <pair <Int_t,Int_t > > Ids; // GlobId, DetId

            while(true) {
                Long64_t TestTS = (Long64_t)testevt->time();

                Long64_t DT = TestTS - CurrentTS;

                if(DT > fDTClover) break;

                if((testdet->IsClover() || testdet->IsAC()) && testdet->CloverId == aclov.ClovID && testevt->is_in_trigger) {

                    if(header.timestamp == 0) header.timestamp = TestTS;

                    if(testdet->IsClover()) {
                        if(DT>0) fDTvsEHists[det->GlobId]->Fill(TestTS-header.timestamp,testevt->data[0]);

                        aclov.Ge_dt[testdet->Id]     = (UShort_t)(TestTS-header.timestamp);
                        aclov.Ge_E[testdet->Id]      = (UShort_t)testevt->data[0];
                        aclov.Ge_Extra[testdet->Id]  = (UShort_t)testevt->data[1];
                        if(fPrintEvts) GeIds += Form(" %d (%d)",testdet->adc,aclov.Ge_dt[testdet->Id]);

                        if(testevt!=CurrentEvt) fDTGeToFirstGe->Fill(DT);
                        if(testdet->IsFIPPS()) fNGoodEvents[1]++;
                        else if(testdet->IsIFIN()) fNGoodEvents[2]++;

                        Ids.push_back(make_pair(testdet->GlobId,testdet->Id));
                    }
                    else if(testdet->IsAC() && DT > fDTACMin && DT<fDTACMax){
                        aclov.AC_dt[testdet->Id]     = (UShort_t)(TestTS-header.timestamp);
                        aclov.AC_E[testdet->Id]      = (UShort_t)testevt->data[0];
                        aclov.AC_Extra[testdet->Id]  = (UShort_t)testevt->data[1];

                        if(fPrintEvts) ACIds += Form(" %d (%d)",testdet->adc,aclov.AC_dt[testdet->Id]);
                        fDTACToFirstGe->Fill(DT);
                        if(testdet->IsFIPPSAC()) fNGoodEvents[5+testdet->Id]++;
                        else if(testdet->IsIFINAC()) fNGoodEvents[3]++;

                        ACMultInClov++;
                    }

                    testevt->is_in_trigger = false;
                }

                testpos++;
                testevt = &(*(testpos));
                if(testpos == last) break;
                testdet = fDetectors.at(testevt->adc);
            };

            fOutputFile.write((char*)&header, sizeof(DetHeader));
            fOutputFile.write((char*)&aclov, sizeof(Clov));

            fNGoodEvents[0]++;

            if(fPrintEvts) {
                cout << left << " ==> WRITE CLOV => Clov Id:" << setw(2) << aclov.ClovID << "First TS: " << header.timestamp << endl;
                cout         << "     Ge adc(dT):" << GeIds << endl;
                cout         << "     AC adc(dT):" << ACIds << endl;
                cout         << "     ClovMult: " << Ids.size() + ACMultInClov  << endl;
            };

            fGeMult->Fill(Ids.size());
            if(ACMultInClov==0)
                for(int i=0 ; i<Ids.size() ; i++) {
                    if(fPlotM1 && Ids.size()==1) fRawGeHists_M1[Ids[i].first]->Fill(aclov.Ge_E[Ids[i].second]);
                    if(fPlotM2 && Ids.size()==2) fRawGeHists_M2[Ids[i].first]->Fill(aclov.Ge_E[Ids[i].second]);
                    if(fPlotM3 && Ids.size()==3) fRawGeHists_M3[Ids[i].first]->Fill(aclov.Ge_E[Ids[i].second]);
                    if(fPlotM4 && Ids.size()==4) fRawGeHists_M4[Ids[i].first]->Fill(aclov.Ge_E[Ids[i].second]);
                }
        }
        else if(det->IsTAG()) {
            DetHeader header;
            TAG atag;
            header.type = kTAG;
            header.timestamp = CurrentTS;
            atag.TAGId = det->CloverId;
            atag.QTot = CurrentEvt->data[0];
            atag.QShort = CurrentEvt->data[1];

            fOutputFile.write((char*)&header, sizeof(DetHeader));
            fOutputFile.write((char*)&atag, sizeof(TAG));

            fNGoodEvents[4]++;

            if(fPrintEvts) {
                cout << left << setw(15) << "New TAG Evt: " << setw(10) << header.timestamp << setw(3) << atag.TAGId  << setw(10) << atag.QTot << setw(10) << atag.QShort << endl;
            }
        }
        else if(det->IsLaBr3()) {
            DetHeader header;
            LaBr3 alabr3;
            header.type = kLaBr3;
            header.timestamp = CurrentTS;
            alabr3.Id = det->CloverId;
            alabr3.energy = CurrentEvt->data[0];
            alabr3.extra = CurrentEvt->data[1];

            fOutputFile.write((char*)&header, sizeof(DetHeader));
            fOutputFile.write((char*)&alabr3, sizeof(LaBr3));

            fNGoodEvents[8]++;

            if(fPrintEvts) {
                cout << left << setw(15) << "New LaBr3 Evt: " << setw(10) << header.timestamp << setw(3) << alabr3.Id  << setw(10) << alabr3.energy << setw(10) << alabr3.extra << endl;
            }
        }
        else if(det->IsTAC()) {
            DetHeader header;
            TAC atac;
            header.type = kTAC;
            header.timestamp = CurrentTS;
            atac.Id = det->CloverId;
            atac.time = CurrentEvt->data[0];
            atac.extra = CurrentEvt->data[1];

            fOutputFile.write((char*)&header, sizeof(DetHeader));
            fOutputFile.write((char*)&atac, sizeof(TAC));

            fNGoodEvents[9]++;

            if(fPrintEvts) {
                cout << left << setw(15) << "New TAC Evt: " << setw(10) << header.timestamp << setw(3) << atac.Id  << setw(10) << atac.time << setw(10) << atac.extra << endl;
            }
        }
        else if(det->IsSiPM()) {
            DetHeader header;
            SiPM sipm;
            header.type = kSiPM;
            header.timestamp = CurrentTS;
            sipm.Id = det->CloverId;
            sipm.q_long = CurrentEvt->data[0];
            sipm.q_short = CurrentEvt->data[1];
            sipm.fine_ts = CurrentEvt->fine_timestamp;

            fOutputFile.write((char*)&header, sizeof(DetHeader));
            fOutputFile.write((char*)&sipm, sizeof(SiPM));

            fNGoodEvents[10]++;

            if(fPrintEvts) {
                cout << left << setw(15) << "New SiPM Evt: " << setw(10) << header.timestamp << setw(3) << sipm.Id  << setw(10) << sipm.q_long << setw(10) << sipm.q_short << setw(10) << sipm.fine_ts << endl;
            }
        }

        ++progress;
        position++;
    }

    cout<<endl;

    return last;
}

void WriteHistos()
{
    // Ana time alignment
    for(int i=1 ; i<= fTSDiffvsIdAlign->GetNbinsY() ; i++) {
        TH1 *proj = fTSDiffvsIdAlign->ProjectionX("proj",i,i);
        Int_t dt = 0;
        if(proj->GetEntries()>100) dt = proj->GetBinCenter(proj->GetMaximumBin());

        if(TMath::Abs(dt)>30) {
            WARN_MESS << "Time alignment for id " << fTSDiffvsIdAlign->GetYaxis()->GetBinLowEdge(i) << " too large: " << dt << "ns" << endl << "Press a key to continue." << ENDL;
            if(fBatchMode == false) cin.get();
        }
        delete proj;
    }

    fOutputFileROOT->cd();

    fOutputFileROOT->mkdir("TimeAlignment");
    fOutputFileROOT->cd("TimeAlignment");
    fTSDiffvsIdRaw->Write();
    fTSDiffvsIdAlign->Write();

    fOutputFileROOT->mkdir("Delta_TS");
    fOutputFileROOT->cd("Delta_TS");
    fDeltaTSvsId->Write();
    fDeltaTS->Write();
    fDeltaTS_Ge->Write();
    fDeltaTS_AC->Write();
    if(fUseTAG) fDeltaTS_TAG->Write();
    if(fUseLaBr3) fDeltaTS_LaBr3->Write();
    if(fUseTAC) fDeltaTS_TAC->Write();
    if(fUseSiPM) fDeltaTS_TSiPM->Write();

    fOutputFileROOT->mkdir("EGe");
    fOutputFileROOT->mkdir("DTvsE");
    if(fPlotM1) fOutputFileROOT->mkdir("EGe_M1");
    if(fPlotM2) fOutputFileROOT->mkdir("EGe_M2");
    if(fPlotM3) fOutputFileROOT->mkdir("EGe_M3");
    if(fPlotM4) fOutputFileROOT->mkdir("EGe_M4");

    fOutputFileROOT->mkdir("AC");

    if(fUseTAG) fOutputFileROOT->mkdir("TAG");
    if(fUseLaBr3) fOutputFileROOT->mkdir("LaBr3");
    if(fUseTAC) fOutputFileROOT->mkdir("TAC");
    if(fUseSiPM) fOutputFileROOT->mkdir("SiPM");

    for (auto const& entry: fDetectorsFromIds) {
        DetDef *det = entry.second;
        if(det->IsClover()) {

            fOutputFileROOT->cd("DTvsE");
            fDTvsEHists[det->GlobId]->Write();

            fOutputFileROOT->cd("EGe");
            fGeHists[det->GlobId]->Write();

            if(fPlotM1){fOutputFileROOT->cd("EGe_M1");fRawGeHists_M1[det->GlobId]->Write();}
            if(fPlotM2){fOutputFileROOT->cd("EGe_M2");fRawGeHists_M2[det->GlobId]->Write();}
            if(fPlotM3){fOutputFileROOT->cd("EGe_M3");fRawGeHists_M3[det->GlobId]->Write();}
            if(fPlotM4){fOutputFileROOT->cd("EGe_M4");fRawGeHists_M4[det->GlobId]->Write();}
        }
        else if(det->IsAC()) {
            fOutputFileROOT->cd("AC");
            fACHists[det->GlobId]->Write();
        }
        else if(fUseTAG && det->IsTAG()) {
            fOutputFileROOT->cd("TAG");
            fTAG_QTot[det->Id]->Write();
            fTAG_QShort[det->Id]->Write();
        }
        else if(fUseLaBr3 && det->IsLaBr3()) {
            fOutputFileROOT->cd("LaBr3");
            fLaBr3_E[det->Id]->Write();
        }
        else if(fUseTAC && det->IsTAC()) {
            fOutputFileROOT->cd("TAC");
            fTAC_T[det->Id]->Write();
        }
        else if(fUseSiPM && det->IsSiPM()) {
            fOutputFileROOT->cd("SiPM");
            fSiPM_QTot[det->Id]->Write();
            fSiPM_QShort[det->Id]->Write();
            fSiPM_T[det->Id]->Write();
        }
    }


    fOutputFileROOT->mkdir("ClovEvts");
    fOutputFileROOT->cd("ClovEvts");
    fGeMult->Write();
    fDTGeToFirstGe->Write();
    fDTACToFirstGe->Write();
    fGeMult = new TH1F("GeMult","GeMult)",5,0,5);

    cout<<endl;
    INFO_MESS << "Output ROOT file written  : " << fOutputFileROOT->GetName() << ENDL ;

    fOutputFileROOT->Close();
}

void InitHistos()
{
    fOutputFileROOT = new TFile(fOutputFileName.Copy().ReplaceAll(".sort",".root"),"recreate");

    fTSDiffvsIdRaw   = new TH2F("TSDiff_vs_Id_Raw","TSDiff_vs_Id_Raw",fTimeAlignmentWindow*2/fTimeBinning,-fTimeAlignmentWindow,fTimeAlignmentWindow,fNIds+1,0,fNIds+1);
    fTSDiffvsIdAlign = new TH2F("TSDiff_vs_Id_Align","TSDiff_vs_Id_Align",fTimeAlignmentWindow*2/fTimeBinning,-fTimeAlignmentWindow,fTimeAlignmentWindow,fNIds+1,0,fNIds+1);

    map < Int_t, DetDef* >::iterator itr;
    for (itr = fDetectorsFromIds.begin(); itr != fDetectorsFromIds.end(); ++itr) {
        DetDef *det = itr->second;
        if(det->IsClover()) {
            if(fPlotM1) fRawGeHists_M1[det->GlobId] = new TH1D(Form("ERawGe_M1_Id%d_ADC%d",det->GlobId,det->adc),Form("ERawGe_M1_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);
            if(fPlotM2) fRawGeHists_M2[det->GlobId] = new TH1D(Form("ERawGe_M2_Id%d_ADC%d",det->GlobId,det->adc),Form("ERawGe_M2_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);
            if(fPlotM3) fRawGeHists_M3[det->GlobId] = new TH1D(Form("ERawGe_M3_Id%d_ADC%d",det->GlobId,det->adc),Form("ERawGe_M3_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);
            if(fPlotM4) fRawGeHists_M4[det->GlobId] = new TH1D(Form("ERawGe_M4_Id%d_ADC%d",det->GlobId,det->adc),Form("ERawGe_M4_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);

            fGeHists[det->GlobId] = new TH1D(Form("Ge_Id%d_ADC%d",det->GlobId,det->adc),Form("Ge_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);

            fDTvsEHists[det->GlobId] = new TH2D(Form("DTvsE%d_ADC%d",det->GlobId,det->adc),Form("DTvsE%d_ADC%d",det->GlobId,det->adc),fTimeAlignmentWindow/fTimeBinning,0,fTimeAlignmentWindow,5000,0,5000);
        }
        else if(det->IsAC()) {
            fACHists[det->GlobId] = new TH1D(Form("AC_Id%d_ADC%d",det->GlobId,det->adc),Form("AC_Id%d_ADC%d",det->GlobId,det->adc),32768,0,32768);

        }
        else if(fUseTAG && det->IsTAG()) {
            fTAG_QTot[det->Id]   = new TH1D(Form("TAG_%d_QTot_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("TAG%d_QTot_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
            fTAG_QShort[det->Id] = new TH1D(Form("TAG_%d_QShort_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("TAG%d_QShort_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
        }
        else if(fUseLaBr3 && det->IsLaBr3()) {
            fLaBr3_E[det->Id]   = new TH1D(Form("LaBr3_%d_E_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("LaBr3%d_E_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
        }
        else if(fUseTAC && det->IsTAC()) {
            fTAC_T[det->Id]   = new TH1D(Form("TAC_%d_T_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("TAC%d_T_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
        }
        else if(fUseSiPM && det->IsSiPM()) {
            fSiPM_T[det->Id]   = new TH1D(Form("SiPM_%d_T_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("SiPM%d_T_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
            fSiPM_QTot[det->Id]   = new TH1D(Form("SiPM_%d_QTot_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("SiPM%d_QTot_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
            fSiPM_QShort[det->Id]   = new TH1D(Form("SiPM_%d_QShort_Id%d_ADC%d",det->Id,det->GlobId,det->adc),Form("SiPM%d_QShort_Id%d_ADC%d",det->Id,det->GlobId,det->adc),32768,0,32768);
        }
    }

    fDeltaTSvsId = new TH2F("DeltaTSvsId","DeltaTSvsId (ns)",1000,0,10000,fNIds+1,0,fNIds+1);
    fDeltaTS     = new TH1F("DeltaTS","DeltaTS (ns)",1000,0,10000);
    fDeltaTS_Ge  = new TH1F("DeltaTS_Ge","DeltaTS_Ge (ns)",1000,0,10000);
    if(fUseTAG) {
        fDeltaTS_TAG = new TH1F("DeltaTS_TAG","DeltaTS_TAG (ns)",1000,0,10000);
    }
    if(fUseLaBr3) {
        fDeltaTS_LaBr3 = new TH1F("DeltaTS_LaBr3","DeltaTS_LaBr3 (ns)",1000,0,10000);
    }
    if(fUseTAC) {
        fDeltaTS_TAC = new TH1F("DeltaTS_TAC","DeltaTS_TAC (ns)",1000,0,10000);
    }
    if(fUseSiPM) {
        fDeltaTS_TSiPM = new TH1F("DeltaTS_SiPM","DeltaTS_SiPM (ns)",1000,0,10000);
    }
    fDeltaTS_AC  = new TH1F("DeltaTS_AC","DeltaTS_AC (ns)",1000,0,10000);

    // Clover evts

    fGeMult = new TH1F("GeMult","GeMult)",10,-0.5,9.5);
    fDTGeToFirstGe = new TH1F("DTGeToFirstGe","DTGeToFirstGe (ns)",fDTClover*2/fTimeBinning,0.,2*fDTClover);
    fDTACToFirstGe = new TH1F("DTACToFirstGe","DTACToFirstGe (ns)",fDTClover*4/fTimeBinning,-2*fDTClover,2*fDTClover);
}
