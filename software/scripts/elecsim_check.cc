#include "Event/ElecHeader.h"
#include "Event/ElecTruthHeader.h"
#include "Context/TimeStamp.h"

void elecsim_check()
{


    //TFile* infile = new TFile("/junofs/users/yumiao/simulation/software_test/trunk/newtrunk/gamma-mix/elecsim-1011.root", "read");
    TFile* infile = new TFile("/junofs/users/yumiao/simulation/software_test/trunk/newtrunk/SN-mix/elecsim-0.root", "read");
    TTree* tElec  = (TTree*)infile->Get("/Event/Elec/ElecEvent");
    TTree* tTruth = (TTree*)infile->Get("/Event/Sim/Truth/LpmtElecTruthEvent");

    JM::ElecEvent* ee = new JM::ElecEvent();
    JM::LpmtElecTruthEvent* et = new JM::LpmtElecTruthEvent();
    tElec->SetBranchAddress("ElecEvent", &ee);
    tTruth->SetBranchAddress("LpmtElecTruthEvent", &et);
    cout << "Total ElecEvents : " << tElec->GetEntries() << endl;
    cout << "Total ElecTruth Events : " << tTruth->GetEntries() << endl;
    if(tElec->GetEntries() != tTruth->GetEntries()) cout << "Wrong input data, please check !" << endl;



    // output definition : 
    std::string pulseTag[5] = { "kNormalPulse", "kDarkPulse", "kAfterPulse", "kDNAfterPulse", "kUnknownPulse" };
    TFile* out = new TFile("checkElecSim-SNmix.root", "recreate");
    TTree* mytree = new TTree("mytree", "elecsim check tree");
    int m_evtID;
    int m_nNormalPulse;
    int m_nDarkPulse;
    int m_nAfterPulse;
    int m_nDNAfterPulse;
    int m_nPulse;
    int m_tag[100000];
    double m_pulsetime[100000];
    double m_pulseamp[100000];
    double m_tts[100000];
    mytree->Branch("evtID", &m_evtID, "evtID/I");
    mytree->Branch("nPulse", &m_nPulse, "nPulse/I");
    mytree->Branch("nNormalPulse", &m_nNormalPulse, "nNormalPulse/I");
    mytree->Branch("nDarkPulse", &m_nDarkPulse, "nDarkPulse/I");
    mytree->Branch("nAfterPulse", &m_nAfterPulse, "nAfterPulse/I");
    mytree->Branch("nDNAfterPulse", &m_nDNAfterPulse, "nDNAfterPulse/I");
    mytree->Branch("hittime", m_pulsetime, "hittime[nPulse]/D");
    mytree->Branch("tag", m_tag, "tag[nPulse]/I");
    mytree->Branch("amplitude", m_pulseamp, "amplitude[nPulse]/D");
    mytree->Branch("TTS", m_tts, "TTS[nNormalPulse]/D");



    //for(int i=0; i<1; i++) {
    for(int i=0; i<tElec->GetEntries(); i++) {
        m_evtID = i ;
        tElec->GetEntry(i);
        tTruth->GetEntry(i);

        // read elecsim pulse truth info
        //const JM::ElecFeeCrate& efc = ee->elecFeeCrate();
        //JM::ElecFeeCrate* m_crate;
        //m_crate = const_cast<JM::ElecFeeCrate*>(&efc);
        //std::map<int, JM::ElecFeeChannel>* feeChannels = &m_crate->channelData();
        //// std::map<int, JM::ElecFeeChannel> &feeChannels = (m_crate->channelData());
        //std::map<int, JM::ElecFeeChannel>::iterator it;


        // read elecsim truth:
        m_nNormalPulse = 0;
        m_nDarkPulse = 0;
        m_nAfterPulse = 0;

        const std::vector<JM::LpmtElecTruth>& let = et->truths();
        std::vector<JM::LpmtElecTruth>* m_truths;
        m_truths = const_cast<std::vector<JM::LpmtElecTruth>*>(&let);
        m_nPulse = m_truths->size();
        //cout << m_truths->size() << endl;
        for (int j = 0; j < m_truths->size(); j++) {
            std::string m_pulsetype = m_truths->at(j).pulsetype();      // sourse type of the pulse('tag' setting  by user, DN, AP or Unknow) 
            int m_pmtId = m_truths->at(j).pmtId();                      // pmt id of the pulse
            int m_npe = m_truths->at(j).npe();                          // photon-electron number of the pulse
            double m_hitTime = m_truths->at(j).hitTime();               // hit time as from detector simulation
            double m_amplitude = m_truths->at(j).amplitude();           // amplitude as from PMT simulation
            double m_TTS = m_truths->at(j).tts();                       // transit time of the hit
            double m_timeoffset = m_truths->at(j).timeoffset();         // time offset of the hit
            TimeStamp m_pulseHitTime = m_truths->at(j).pulseHitTime();  // TimeStamp of the pulse = event time stamp + transit time + time offset
            if(m_pulsetype == "SN" ) {
                m_tag[j] = 0;
                m_tts[m_nNormalPulse] = m_TTS;
                m_nNormalPulse++;
            } else if (m_pulsetype == "DarkPulse" ) {
                m_nDarkPulse++;
                m_tag[j] = 1;
            } else if(m_pulsetype == "AfterPulse") {
                m_nAfterPulse++;
                m_tag[j] = 2;
            }
            m_pulsetime[j] = m_hitTime;
            m_pulseamp[j]  = m_amplitude;
        
        }

        out->cd();
        mytree->Fill();

    }

    mytree->Write();
    out->Close();
}
