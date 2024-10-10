#include <iostream>
#include <iomanip>

#include "TRint.h"
#include "TGClient.h"
#include "TString.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TSystem.h"

using namespace std;

TRint *theApp;

int main(int argc, char **argv)
{
    if(getenv("FIPPS_Softs_SYS") == nullptr) {
        cout<<"Environment variable: FIPPS_Softs_SYS needs to be define ==> EXIT"<<endl;
        return 1;
    }

    gSystem->Load("libMathMore.so");

    gSystem->CompileMacro(Form("%s/include/AngCorrLib.h",getenv("FIPPS_Softs_SYS")),"dg");
    gSystem->AddIncludePath(Form("-I %s/include",getenv("FIPPS_Softs_SYS")));

    gEnv->SetValue("Gui.IconPath",Form("%s/icons:%s/icons",getenv("ROOTSYS"),getenv("GWSYS")));
    gEnv->SetValue("X11.UseXft","yes");
    gEnv->SetValue("Canvas.ShowGuideLines","false");
    gEnv->SetValue("Unix.*.Root.UseTTFonts","true");

    cout << "  **************************************************************\n";
    cout << "  *    HELLO  -- You are Running Angular Correlation program   *\n";
    cout << "  **************************************************************\n";

    TString file="";
    if(argc==2) {
        file = argv[1];
    }

    theApp = new TRint("App", &argc, argv,nullptr,0);

    // make sure that the Gpad and GUI libs are loaded
    TRint::NeedGraphicsLibs();
    TApplication::NeedGraphicsLibs();
    gApplication->InitializeGraphics();

    theApp->Run(false);

    if (theApp) {
        delete (theApp);
        theApp = nullptr;
    } 
    return 0;
}
