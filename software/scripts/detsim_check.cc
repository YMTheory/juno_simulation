#include "Event/SimHeader.h"
#include  "Event/SimEvent.h"
#include "Event/SimTrack.h"

void detsim_check()
{
    TFile *file = new TFile("/junofs/users/yumiao/simulation/software_test/J20v2r0-branch/detsim-6002.root", "read");
    TTree *t_simEvent = (TTree*)file->Get("/Event/Sim/SimEvent");
    
    std::vector<JM::SimTrack*> vec_simTack;
    t_simEvent->SetBranchAddress("m_tracks", &vec_simTack);
    Int_t nhits;
    t_simEvent->SetBranchAddress("m_nhits", &nhits);
    Int_t ntrks;
    t_simEvent->SetBranchAddress("m_ntrks",&ntrks);

    cout << "Total entry in file: " << t_simEvent->GetEntries() << endl;
    for(int i=0; i<10/*t_simEvent.GetEntries()*/; i++) {
        t_simEvent->GetEntry(i);
        cout << "nTracks : " << nhits << endl;
        for(int itrack=0; itrack<vec_simTack.size(); itrack++) {
            cout << vec_simTack[itrack]->getPDGID();
        }
    }


    delete t_simEvent;
    delete file;
}
