#ifndef PMTSimAlg_h
#define PMTSimAlg_h
#include "SniperKernel/AlgBase.h"
#include "ElecSimSvc.h"
#include "IPMTTool.h"
#include "EvtNavigator/NavBuffer.h"
#include "InputReviser/InputReviser.h"
#include "Event/SimHeader.h"
#include "EventKeeper.h"
#include "SpmtElecConfigSvc.h"
#include "SniperKernel/Task.h" 


using namespace std;
#define NUM_BINS_DIST 96
static double afterPusleTimingDist[NUM_BINS_DIST] = {0.0000000, 0.0012950, 0.0030194, 0.0064297, 0.0142225,
    0.0247712, 0.0417949, 0.0696709, 0.1063259, 0.1485606,
    0.1957973, 0.2465997, 0.2957179, 0.3401377, 0.3798086,
    0.4135611, 0.4437462, 0.4702155, 0.4928434, 0.5123140,
    0.5289294, 0.5432514, 0.5556360, 0.5664919, 0.5758642,
    0.5843533, 0.5920903, 0.5990121, 0.6052796, 0.6111246,
    0.6165560, 0.6215807, 0.6261028, 0.6304032, 0.6346304,
    0.6387545, 0.6425872, 0.6465577, 0.6500606, 0.6534989,
    0.6568080, 0.6602148, 0.6635239, 0.6667457, 0.6700025,
    0.6737148, 0.6780693, 0.6829002, 0.6879599, 0.6935379,
    0.6994667, 0.7058650, 0.7121394, 0.7187122, 0.7255085,
    0.7324600, 0.7394151, 0.7467314, 0.7546551, 0.7630884,
    0.7718917, 0.7816427, 0.7914601, 0.8023700, 0.8137512,
    0.8255861, 0.8375781, 0.8498354, 0.8621276, 0.8741685,
    0.8855410, 0.8962833, 0.9061268, 0.9156126, 0.9243357,
    0.9320726, 0.9392022, 0.9459897, 0.9519587, 0.9574180,
    0.9623468, 0.9667851, 0.9708308, 0.9744296, 0.9777980,
    0.9807231, 0.9835296, 0.9861040, 0.9885404, 0.9906504,
    0.9926209, 0.9943923, 0.9959998, 0.9974658, 0.9987678,
    1.0000000};


class PMTSimAlg: public AlgBase {
 public:
    PMTSimAlg(const std::string& name);
    ~PMTSimAlg();

    bool initialize();
    bool execute();
    bool finalize();

 private:
    int generate_cd_pulses(JM::SimEvent* evt,TimeStamp t);
    int generate_wp_pulses(JM::SimEvent* evt,TimeStamp t);
    int generate_tt_pulses(JM::SimEvent* evt,TimeStamp t);

    int generate_cd_dark_pulses(TimeStamp t1,TimeStamp t2);
    int merge_cd_LPMT_pulses();
    int generate_wp_dark_pulses(TimeStamp t1,TimeStamp t2);
    int merge_wp_LPMT_pulses();
 
 private:
    // cd Lpmt relative
    double PulseAmp(double weight,double gain, double sigmaGain);
    double PulseAmpMCP(double weight,double gain, double sigmaGain, double exp_frac);
    double ConvertPdfRand01(double rand,std::vector<double> pdf, std::vector<double> edges);
    int PoissonRand(double mean);

 private:
    std::shared_ptr<JM::EvtNavigator> get_one_event();
    std::shared_ptr<JM::EvtNavigator> get_one_SN_event();
    EventKeeper::Entry m_current_event;
    deque<EventKeeper::Entry> vec_current_event;
 private:
    ElecSimSvc* m_elecsvc;
    std::map<std::string, Task*> m_taskObjMap;

    SpmtElecConfigSvc* m_SpmtConfigSvc;

    std::deque<std::shared_ptr<JM::EvtNavigator>> vecNav;
    std::deque<TimeStamp> vecEvtTime;
    std::deque<double> vecTimeWindow;

    // if PMT categories are quite different, we just use
    // different PMT tools to generate pulse.
    std::map<PMTType, IPMTTool*>  m_pmt_tool;

    // load events from different tasks.
    std::map<std::string, std::string> m_taskMap;
    std::map<std::string, InputReviser*> m_incidentMap;

    double m_totalRate;///the sum of events rate
    double m_mainTau;///1.0/totalRate
    double m_darkNoiseWindows;
    bool m_loop;
    std::map<std::string,double> m_rateMap;///events rate map
    std::map<std::string,int> m_startIndexMap;///events start index map
    std::map<std::string,int> m_loopMap;//Event circulation map,first is sample name, second is 0 or 1.

    //used to decide first entry's position
    std::map<std::string, bool> m_firstMap;

    // time related
    TimeStamp curEvtTimeStamp;
    TimeStamp curEvtStartHitTime;
    TimeStamp curEvtEndHitTime;

    TimeStamp nextEvtTime;

    // used by SN
    TimeStamp nextBGKEvtTime;
    TimeStamp nextSNEvtTime;
    TimeStamp firstSNEvtTime;

    bool isSNEvt;

    // pmt numbers
    int m_total_Lpmt;
    int m_total_Spmt;
    int m_total_WPLpmt;
    int m_init_seed;   
    // LPMT option
    bool m_enableLPMT;
    bool m_enableSPMT;
    bool m_enableWP;
    bool m_enableSNMode;
    bool m_enableEfficiency; 
    bool m_enableAfterPulse;
    bool m_enableDarkPulse;
    bool m_enablePmtTTS;
    bool m_enableWPDarkPulse;
    bool m_enableMergeLPMTPulse;
    bool m_enableAssignGain;
    bool m_enableAssignSigmaGain;
    bool m_enableAssignDarkPulse;
    bool m_enableMCPLargeSignal; 
    double m_Gain;
    double m_SigmaGain;
    double m_darkRate;
    double m_darkRate_WP;
    double m_darkRateScale;
    double m_darkRateScale_WP;
    double m_darkRateScale_SPMT;
    float m_evtBufferLength;
    float m_nHitsThreshold;
    std::vector<double> m_afterPulsePdf;
    std::vector<double> m_afterPulseEdges;
    float edepX;
    float edepY;
    float edepZ;

    string m_SN_sampleName;
 
    bool m_stoppingRun;

    TGraph* gNNVTAmp;

    // In order to stop the simulation without loss events, an error counter is added here.
    // if the times failed to load is more than N_LIMIT_COUNTS, stop.
    // -- Tao Lin <lintao@ihep.ac.cn>, 2020.05.26
    int m_count_errs;
    int m_nlimit_counts;

    //add for SN Mode sim
    std::shared_ptr<JM::EvtNavigator> SN_aNav;
    std::deque<std::shared_ptr<JM::EvtNavigator>> bufferNav;
};

#endif
