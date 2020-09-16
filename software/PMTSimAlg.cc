#include "PMTSimAlg.h"
#include <boost/python.hpp> 
#include <ios> 

#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolBase.h"
#include "SniperKernel/Incident.h"
#include "SniperKernel/SniperException.h"

#include "InputReviser/InputReviser.h"
#include "SniperKernel/Task.h"

#include "Event/SimHeader.h"

#include <TTimeStamp.h>
#include "Context/TimeStamp.h"

#include <TRandom.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include "SPMTTool.h"
#include "TF1.h"


DECLARE_ALGORITHM(PMTSimAlg);

using namespace std;

PMTSimAlg::PMTSimAlg(const std::string& name)
    : AlgBase(name), m_elecsvc(0)
{
    declProp("Seed", m_init_seed);
    declProp("RateMap", m_rateMap);
    declProp("TaskMap", m_taskMap);
    declProp("LoopMap", m_loopMap);
    declProp("StartIdxMap", m_startIndexMap);


    m_total_Lpmt = 17613;
    m_total_Spmt = 25600;
    m_total_WPLpmt =2400;

    declProp("TotalLpmt",m_total_Lpmt=17613);
    declProp("enableLPMT", m_enableLPMT=true);
    declProp("enableSPMT", m_enableSPMT=true);
    declProp("enableWP", m_enableWP=false);
    declProp("enableSNMode", m_enableSNMode=false);
    declProp("enableEfficiency", m_enableEfficiency=false);
    declProp("enableAfterPulse", m_enableAfterPulse=false);
    declProp("enableDarkPulse", m_enableDarkPulse=false);
    declProp("enablePmtTTS", m_enablePmtTTS=true);
    declProp("enableWPDarkPulse", m_enableWPDarkPulse=false);
    declProp("enableMergeLPMTPulse", m_enableMergeLPMTPulse=false);
    declProp("enableAssignGain", m_enableAssignGain=false);
    declProp("enableAssignDarkPulse", m_enableAssignDarkPulse=false);
    declProp("enableAssignSigmaGain", m_enableAssignSigmaGain=false);
    declProp("enableMCPLargeSignal", m_enableMCPLargeSignal=true);
    declProp("inputGain", m_Gain=1);
    declProp("inputSigmaGain", m_SigmaGain=0.3);
    declProp("darkRate", m_darkRate=30e3);
    declProp("darkRateScale", m_darkRateScale=1.);
    declProp("darkNoiseWindows", m_darkNoiseWindows=15e5);
    declProp("darkRate_WP", m_darkRate_WP=70e3);
    declProp("darkRateScale_WP", m_darkRateScale_WP=1.);
    declProp("darkRateScale_SPMT", m_darkRateScale_SPMT=1.);
    declProp("evtBufferLength", m_evtBufferLength=5000);
    declProp("nHitsThreshold", m_nHitsThreshold=200);

    m_afterPulsePdf.clear();
    m_afterPulseEdges.clear();

    for(int ii=0; ii < NUM_BINS_DIST; ii++) {
        m_afterPulsePdf.push_back(afterPusleTimingDist[ii]);
        m_afterPulseEdges.push_back(ii*100.+500.);
    }

    declProp("N_LIMIT_COUNTS", m_nlimit_counts=3);
    m_count_errs = 0;

    m_stoppingRun = false;
}

PMTSimAlg::~PMTSimAlg()
{

}

bool
PMTSimAlg::initialize()
{
    LogInfo<<"Initialize the SEED = "
           <<m_init_seed
           <<std::endl;
    gRandom->SetSeed(m_init_seed);
 
    SniperPtr<ElecSimSvc> svc(*getRoot(), "ElecSimSvc");
    if (svc.invalid()) {
        LogError << "can't find service ElecSimSvc" << std::endl;
        return false;
    }
    m_elecsvc = svc.data();
    curEvtTimeStamp = m_elecsvc->get_start_time();
    nextEvtTime = m_elecsvc->get_start_time();
    nextSNEvtTime = m_elecsvc->get_start_time();
    nextBGKEvtTime = m_elecsvc->get_start_time();

    isSNEvt = true;
    SN_aNav = nullptr;

    gNNVTAmp=m_elecsvc->get_graph();

    // register pmt tool
    // here, we use the same pmt tool, with different configurations.
    m_pmt_tool[kPMT_CD_Lpmt] = tool<IPMTTool>("LPMTTool"); // -> PMTTool
    m_pmt_tool[kPMT_CD_Spmt] = tool<IPMTTool>("SPMTTool"); // -> PMTTool
    m_pmt_tool[kPMT_WP] =      tool<IPMTTool>("LPMTTool"); // -> PMTTool
    // otherwise, we can use the same tool instead
    // m_pmt_tool[kPMT_CD_Lpmt] = tool<IPMTTool>("PMTTool"); // -> PMTTool
    // m_pmt_tool[kPMT_CD_Spmt] = tool<IPMTTool>("PMTTool"); // -> PMTTool

    // check whether hander exists
    for (std::map<PMTType, IPMTTool*>::iterator it = m_pmt_tool.begin();
         it != m_pmt_tool.end(); ++it) {
        if (!it->second) {
            LogError << "Can't find pmt tool for " << it->first << std::endl;
            return false;
        }
    }

    m_SN_sampleName = "";

    m_totalRate = 0;
    for(map<string,double>::iterator it=m_rateMap.begin();it!=m_rateMap.end();it++){ //key is sample name ,value is the rate of sample
        m_firstMap.insert(std::pair<string,bool>(it->first,true));// I define this map just for find the begin entry of the sample
        //check weather SN is enabled, if yes, reset the rate of SN to be zero.
        if(m_enableSNMode) {
          if(strstr(it->first.c_str(), "SN") == NULL) m_totalRate += it->second;  //first is sample name , second is rate 
          else {
            m_SN_sampleName = it->first;
            it->second = 0;
          }
        }
        else m_totalRate+=it->second;  //first is sample name , second is rate 
    }
    LogInfo<<"event totalRate= "<< m_totalRate << endl;
    m_mainTau=1.0/m_totalRate;
    LogInfo<< "mainTau= " << m_mainTau <<endl;

    nextBGKEvtTime.Add(gRandom -> Exp(m_mainTau));
                  
    // prepare input reviser
    // m_taskMap["default"] = "default";
    //        ^- sample name      ^- task name

    Task* toptask = getRoot();

    for(std::map<std::string,std::string>::iterator it=m_taskMap.begin();
        it!=m_taskMap.end();it++){
        m_loop = false; // false: don't reuse sample.

        map<string,int>::iterator i=m_loopMap.find(it->first);
        if (i!=m_loopMap.end()) {
            if (i->second==1) m_loop = true;else if(i->second==0) m_loop = false;else m_loop = true;
        }
        Task* task = dynamic_cast<Task*>(toptask->find(it->second));
        if (!task) {
            LogError << "can't find task " << it->second << std::endl;
            return false;
        }
        InputReviser* aInci = new InputReviser(*toptask, it->second, m_loop);  //first is task name 
        m_taskObjMap[it->first] = task;
        m_incidentMap.insert(make_pair(it->first, aInci));//first is sample name，second is the InputReviser which related to the sample name. In test I use sample as task name.
    }

    // interface to Spmt Configuration
    SniperPtr<SpmtElecConfigSvc> spmtSvc(*getRoot(), "SpmtElecConfigSvc");
    if (spmtSvc.invalid()) {
        LogError << "can't find service SpmtElecConfigSvc" << std::endl;
        return false;
    }
    m_SpmtConfigSvc = spmtSvc.data();
 
    SPMTTool* spmtTool = (SPMTTool *)m_pmt_tool[kPMT_CD_Spmt]; 
    spmtTool->setSpmtConfig(m_SpmtConfigSvc);

    return true;
}

bool
PMTSimAlg::execute()
{
    // event mixing is also here.
    // we will mix the events, then convert the hit into pulse.
    // so we can avoid duplication of SimHit 
    // 

    LogDebug << "In PMTSimAlg, going to load events and unpacking hits." << std::endl;

    std::deque<Pulse>& CDSPMT_pulse_buffer = m_elecsvc->get_CDSPMT_pulse_buffer();
    std::deque<Pulse>& CDLPMT_pulse_buffer = m_elecsvc->get_CDLPMT_pulse_buffer();
    std::deque<Pulse>& WPLPMT_pulse_buffer = m_elecsvc->get_WPLPMT_pulse_buffer();

    int cnt_loaded_evt = 0;

    while(true) {
        //events pre select optimizate simulation speed.
        int hitsSum = 0;
        TimeStamp timeStart, timeCurrent, timeStop, TempT;
        double delta = 0;

        vec_current_event.clear();
        
        //start time of this round of event sampling
        timeStart = nextEvtTime; 

        //event time of the current event
        timeCurrent = nextEvtTime;
     
        LogDebug << "The time stamp of the first event that might be loaded: " << nextEvtTime << std::endl;

        //set the marker time which will be used to indicate the end point of the pulses in buffer to be processed.
        //The events after the marker time will be stopped to be loaded.
        //It will be used by the followed modules, such as TriggerSim.
        TimeStamp markerT;
        if(m_elecsvc->getFirstPulseTimeInBuffer() >= timeCurrent && m_stoppingRun == false) {
          markerT = timeCurrent;
          markerT.Add(m_evtBufferLength*1e-9);
        }
        else {
          markerT = m_elecsvc->getFirstPulseTimeInBuffer();
          markerT.Add(m_evtBufferLength*1e-9);
        }
        m_elecsvc->set_marker_time(markerT);

        LogDebug << "The marker time stamp is: " << markerT << std::endl;
        LogDebug << "First pulse time in buffer: " << m_elecsvc->getFirstPulseTimeInBuffer()  << std::endl;

        if(m_stoppingRun) {
          nextEvtTime = m_elecsvc->get_marker_time();
          nextEvtTime.Add(m_evtBufferLength*1e-9);
        }

        //number of hits in the interested time window, within marker time.
        int nhitsInWoI = 0;

        if(m_elecsvc->getFirstPulseTimeInBuffer() < m_elecsvc->get_marker_time()) {
          for(int ipulse = 0; ipulse < CDLPMT_pulse_buffer.size(); ipulse++) {
            if(CDLPMT_pulse_buffer[ipulse].pulseHitTime < m_elecsvc->get_marker_time()) nhitsInWoI++;
            else break;
          }
        }
        
        LogDebug << "Number of hits in WoI: " << nhitsInWoI << std::endl;
        LogDebug << "Number of hits in loaded event: " << hitsSum << std::endl;

        while( (hitsSum + nhitsInWoI) < m_nHitsThreshold or 
               nextEvtTime < m_elecsvc->get_marker_time()) {
            // skip the low energy events with no overlap with the next sampled event
          if(nextEvtTime > m_elecsvc->get_marker_time() and 
             (hitsSum + nhitsInWoI) < m_nHitsThreshold) {
              while(true) {
                if(CDLPMT_pulse_buffer.size() == 0) {
                  timeStart = nextEvtTime;
                  m_elecsvc->setFirstPulseTimeInBuffer(nextEvtTime);
                  break;
                }
                if(CDLPMT_pulse_buffer[0].pulseHitTime > m_elecsvc->get_marker_time()) {
                  if(CDLPMT_pulse_buffer[0].pulseHitTime > nextEvtTime) timeStart = nextEvtTime;
                  else timeStart = CDLPMT_pulse_buffer[0].pulseHitTime;
                  m_elecsvc->setFirstPulseTimeInBuffer(CDLPMT_pulse_buffer[0].pulseHitTime);
                  break;
                }
                CDLPMT_pulse_buffer.pop_front();
              }
              markerT = timeStart;
              markerT.Add(m_evtBufferLength*1e-9);
              m_elecsvc->set_marker_time(markerT);

                  nhitsInWoI = 0;
                  if(m_elecsvc->getFirstPulseTimeInBuffer() < m_elecsvc->get_marker_time()) {
                    for(int ipulse = 0; ipulse < CDLPMT_pulse_buffer.size(); ipulse++) {
                      if(CDLPMT_pulse_buffer[ipulse].pulseHitTime < m_elecsvc->get_marker_time()) nhitsInWoI++;
                      else break;
                    }
                  }

                  LogDebug << "Marker time is updated to: " << markerT << std::endl;
                  LogDebug << "First pulse time in buffer: " << m_elecsvc->getFirstPulseTimeInBuffer() << std::endl;
                  LogDebug << "Number of pulse in WoI: " << nhitsInWoI << std::endl;

                  vecNav.clear();
                  vecEvtTime.clear();
                  vecTimeWindow.clear();
                  vec_current_event.clear();
                  hitsSum = 0;
              }

              if(nextEvtTime > m_elecsvc->get_marker_time()) continue;

              //the next event is loaded
              std::shared_ptr<JM::EvtNavigator> aNav = nullptr;
              if(m_enableSNMode) {
                if(nextSNEvtTime == m_elecsvc->get_start_time()) { 
                  SN_aNav = get_one_SN_event();
                  if (SN_aNav == nullptr) {
                      LogWarn << "Get empty navigator of SN events. Run is stopping." << std::endl;
                      m_stoppingRun = true;
                      break;
                  }
                  JM::SimHeader* first_hdr = dynamic_cast<JM::SimHeader*>(SN_aNav->getHeader("/Event/Sim"));
                  JM::SimEvent* first_evt = dynamic_cast<JM::SimEvent*>(first_hdr->event());
                  const std::vector<JM::SimTrack*>& trc = first_evt->getTracksVec();
                  std::vector<JM::SimTrack*>::const_iterator it_track = trc.begin();
                  JM::SimTrack* track = *it_track;
                  double track_initT = track->getInitT();

                  TimeStamp evtTime(track_initT*1e-9);
                  firstSNEvtTime = evtTime;
                  aNav = SN_aNav;
                }
                else {
                  if(isSNEvt) aNav = SN_aNav;
                  else aNav = get_one_event();
                }
              }
              else {
                aNav = get_one_event();
              }

              if (aNav == nullptr) {
                  LogWarn << "Get empty navigator. Run is stopping." << std::endl;
                  m_stoppingRun = true;
                  break;
              }

              ++cnt_loaded_evt;

              //to get number of hits of LPMT via simHeader
              JM::SimHeader* hdr = dynamic_cast<JM::SimHeader*>(aNav->getHeader("/Event/Sim"));
              int nHits = hdr->getCDLPMTtotalHits();

              timeCurrent = nextEvtTime;
              aNav->setTimeStamp(TTimeStamp(timeCurrent.GetSec(), timeCurrent.GetNanoSec()));
              vecNav.push_back(aNav); 
              vecEvtTime.push_back(timeCurrent); 
              double tempTimeWindow = hdr->getCDLPMTtimeWindow();
              vecTimeWindow.push_back(tempTimeWindow*1e-9);

              TimeStamp stopTimeOfNextEvt = timeCurrent;
              stopTimeOfNextEvt.Add(tempTimeWindow*1e-9);

              hitsSum += nHits;
              
              LogDebug << "one event is loaded." << std::endl;
              LogDebug << "Time stamp of this event: " << timeCurrent << std::endl;
              LogDebug << "Number of hits in this event: " << hitsSum << std::endl;
              LogDebug << "Number of hits of loaded events: " << hitsSum << std::endl;

              if(m_enableSNMode) {
                if(isSNEvt) {
                  SN_aNav = get_one_SN_event();
                  if (SN_aNav == nullptr) {
                      LogWarn << "Get empty navigator of SN events. Run is stopping." << std::endl;
                      m_stoppingRun = true;
                      break;
                  }
                  JM::SimHeader* first_hdr = dynamic_cast<JM::SimHeader*>(SN_aNav->getHeader("/Event/Sim"));
                  JM::SimEvent* first_evt = dynamic_cast<JM::SimEvent*>(first_hdr->event());
                  const std::vector<JM::SimTrack*>& trc = first_evt->getTracksVec();
                  std::vector<JM::SimTrack*>::const_iterator it_track = trc.begin();
                  JM::SimTrack* track = *it_track;
                  double track_initT = track->getInitT();

                  TimeStamp evtTime(track_initT*1e-9);
                  nextSNEvtTime = m_elecsvc->get_start_time();
                  nextSNEvtTime.Add((evtTime - firstSNEvtTime).GetSeconds());
                }
                else {
                  delta = gRandom -> Exp(m_mainTau);
                  nextBGKEvtTime.Add(delta);
                }

                if(nextBGKEvtTime > nextSNEvtTime) {
                  nextEvtTime = nextSNEvtTime;
                  isSNEvt = true;
                }
                else {
                  nextEvtTime = nextBGKEvtTime;
                  isSNEvt = false;
                }
              }
              else {
                delta = gRandom -> Exp(m_mainTau);
                nextEvtTime = timeCurrent;
                nextEvtTime.Add(delta);
              }
        }

        if(vecEvtTime.size() == 0 && CDLPMT_pulse_buffer.size() == 0) {
            LogWarn << "nothing is loaded." << std::endl;
            break;
        } 
        if(vecEvtTime.size() == 0 && CDLPMT_pulse_buffer.size() > 0) {
            m_elecsvc->set_cur_evttime(m_elecsvc->getFirstPulseTimeInBuffer());
        }
        if(vecEvtTime.size() > 0) {
            m_elecsvc->set_cur_evttime(vecEvtTime.back());
        }

        if(m_enableDarkPulse){
           if(m_elecsvc->getFirstPulseTimeInBuffer() > timeStart) generate_cd_dark_pulses(timeStart, m_elecsvc->get_marker_time());
           else generate_cd_dark_pulses(m_elecsvc->getFirstPulseTimeInBuffer(), m_elecsvc->get_marker_time());
        }

        if(m_enableWPDarkPulse){
          if(m_elecsvc->getFirstPulseTimeInBuffer() > timeStart) generate_wp_dark_pulses(timeStart, m_elecsvc->get_marker_time());
          else generate_wp_dark_pulses(m_elecsvc->getFirstPulseTimeInBuffer(), m_elecsvc->get_marker_time());
        }

        //loop EvtNavigator buffer and genetate pulses
        JM::SimHeader* curr_hdr = 0;
        JM::SimEvent* evt = 0;
        int eventNum = 0;
        for(std::deque<std::shared_ptr<JM::EvtNavigator>>::iterator it_nav = vecNav.begin();it_nav != vecNav.end();++it_nav) {
            std::shared_ptr<JM::EvtNavigator> currNav = *it_nav;
            curr_hdr = dynamic_cast<JM::SimHeader*>(currNav->getHeader("/Event/Sim"));
            evt = dynamic_cast<JM::SimEvent*>(curr_hdr->event());
            //vertex fitting
            JM::SimTrack* tra = 0;
            edepX =0;
            edepY =0;
            edepZ =0;
            const std::vector<JM::SimTrack*>&  tracks = evt->getTracksVec();
            for( std::vector<JM::SimTrack*>::const_iterator it_tra = tracks.begin();
                 it_tra != tracks.end();++it_tra) {
                    tra = *it_tra;
                    edepX = tra->getEdepX();
                    edepY = tra->getEdepY();
                    edepZ = tra->getEdepZ();
            }

            if(abs(edepX)>17500) edepX = edepX/abs(edepX)*15000;
            else edepX = 5000*round(edepX/5000);
            if(abs(edepY)>17500) edepY = edepY/abs(edepY)*15000;
            else edepY = 5000*round(edepY/5000);
            if(abs(edepZ)>17500) edepZ = edepZ/abs(edepZ)*15000;
            else edepZ = 5000*round(edepZ/5000);

            //Event keeper

                   
            EventKeeper& m_keeper = EventKeeper::Instance();
            m_keeper.set_current_entry(vec_current_event[eventNum]);
            
            //generate cd pulses
            generate_cd_pulses(evt,vecEvtTime[eventNum]);  

            if(m_enableWP) generate_wp_pulses(evt,vecEvtTime[eventNum]);

            generate_tt_pulses(evt,vecEvtTime[eventNum]);
            eventNum++;
        }
       
        vecNav.clear();
        vecEvtTime.clear();
        vecTimeWindow.clear();

        std::sort(CDSPMT_pulse_buffer.begin(), CDSPMT_pulse_buffer.end());
        std::sort(CDLPMT_pulse_buffer.begin(), CDLPMT_pulse_buffer.end());

        if(m_enableMergeLPMTPulse){
            merge_cd_LPMT_pulses();
        }
        // need to sort pulse
        std::sort(WPLPMT_pulse_buffer.begin(), WPLPMT_pulse_buffer.end());
        // merge pulses for wp LPMT
        if(m_enableMergeLPMTPulse){
            merge_wp_LPMT_pulses();
        }
        break;
    }

    // if in this event, nothing is loaded,
    // we need to throw the exception again.
    if (!cnt_loaded_evt && CDLPMT_pulse_buffer.size() == 0) {
        LogWarn << "Nothing is loaded, propose to stop the simulation. " << std::endl;
        LogWarn << "The current status of top task in Sniper is "
                << static_cast<unsigned short>(getRoot()->Snoopy().state())
                << std::endl;
        for (const auto& subtask: m_taskObjMap) {
            LogWarn << "The current status of sub task ["
                    << subtask.first << "]: "
                    << static_cast<unsigned short>(subtask.second->Snoopy().state())
                    << std::endl;
        }

        LogWarn << "Number of pulses in CDSPMT: "
                << CDSPMT_pulse_buffer.size()
                << std::endl;
        LogWarn << "Number of pulses in CDLPMT: "
                << CDLPMT_pulse_buffer.size()
                << std::endl;
        LogWarn << "Number of pulses in WPLPMT: "
                << WPLPMT_pulse_buffer.size()
                << std::endl;

        ++m_count_errs;

        if (m_count_errs>=m_nlimit_counts) {
            // mark current task as stop
            getParent()->stop(Sniper::StopRun::Peacefully);
        }
    }

    return true;
}

bool
PMTSimAlg::finalize()
{
    delete gNNVTAmp;
    return true;
}

int 
PMTSimAlg::generate_cd_pulses(JM::SimEvent* evt, TimeStamp curTimeStamp)
{
    std::deque<Pulse>& SPMT_pulse_buffer = m_elecsvc->get_CDSPMT_pulse_buffer();
    std::deque<Pulse>& LPMT_pulse_buffer = m_elecsvc->get_CDLPMT_pulse_buffer();

    int pulses_sum = 0; // generated pulses

    const std::vector<JM::SimPMTHit*>& cd_hits = evt->getCDHitsVec();

    EventKeeper& keeper = EventKeeper::Instance();
    const EventKeeper::Entry& entry = keeper.current_entry();

    for (std::vector<JM::SimPMTHit*>::const_iterator it_hit = cd_hits.begin();
         it_hit != cd_hits.end(); ++it_hit) {

        JM::SimPMTHit* hit=*it_hit;
        int hit_pmtid = hit->getPMTID();
        TimeStamp hitTimeStamp = curTimeStamp;
        TimeStamp corrhitTimeStamp = curTimeStamp;
        double hit_amplitude = 0;
        
        double m_efficiency = 0;
        double m_afterPulseProb = 0;
        double m_timeSpread = 0;
        double m_timeOffset = 0;
        double m_gain = 0;
        double m_sigmaGain = 0;
        double hit_timeSpread = 0;
        double m_tof = 0;
        if(hit_pmtid<=17739){//LPMT
            if(!m_enableLPMT) continue;
            //LPMT data service
            if(m_enableAssignGain){
                m_gain = m_Gain;
            }else{
                m_gain = m_elecsvc->get_gain(hit_pmtid);
            }
            if(m_enableAssignSigmaGain){
                m_sigmaGain = m_SigmaGain;
            }else{
                m_sigmaGain = m_elecsvc->get_sigmaGain(hit_pmtid);
            }

            m_efficiency = m_elecsvc->get_efficiency(hit_pmtid);
            m_afterPulseProb = m_elecsvc->get_afterPulseProb(hit_pmtid);
            m_timeSpread = m_elecsvc->get_timeSpread(hit_pmtid);
            m_timeOffset = m_elecsvc->get_timeOffset(hit_pmtid);

            // here, the hit time is not the absolute timestamp

            if(m_enablePmtTTS) hit_timeSpread =  gRandom->Gaus(0,1) * (m_timeSpread/2.355);
            m_tof = sqrt(pow(edepX-m_elecsvc->get_pmtPosX(hit_pmtid),2)+pow(edepY-m_elecsvc->get_pmtPosY(hit_pmtid),2)+pow(edepZ-m_elecsvc->get_pmtPosZ(hit_pmtid),2))/(2.99792458*1.5*100); 
            double hit_hitTime = hit->getHitTime() + hit_timeSpread + m_timeOffset;

            hitTimeStamp.Add(hit_hitTime*1e-9); // ns -> s  TTS + timeoffset
            double corrtime = hit_hitTime-m_tof;
            if(corrtime<0) corrtime=0;
            corrhitTimeStamp.Add((corrtime)*1e-9);

            double hit_weight = hit->getNPE();
            // add QE affect.
            if(m_enableEfficiency){
                if(gRandom->Rndm() > m_efficiency){
                    continue;    //ignore SimHit due to efficiency
                }
            }

            // single PE resolution
            if(m_elecsvc->get_isHamamatsu(hit_pmtid)==1){          
                hit_amplitude =  PulseAmp(hit_weight,m_gain,m_sigmaGain);
            }else{
                //TGraph* g1 = m_elecsvc->get_graph();
                double theta = hit->getLocalTheta()*180/3.1415926;
                double exp_frac = gNNVTAmp->Eval(theta);
                hit_amplitude =  PulseAmpMCP(hit_weight,m_gain,m_sigmaGain,exp_frac);
            }
        }
        else { //SPMT
            if(!m_enableSPMT) continue;
            double hit_hitTime = hit->getHitTime();
            hitTimeStamp.Add(hit_hitTime*1e-9); // ns -> s  TTS + timeoffset
            hit_amplitude = hit->getNPE();
        }

        IPMTTool* tool = 0;
        if (hit_pmtid<18000) { // LPMT
            tool = m_pmt_tool[kPMT_CD_Lpmt];
        } else if (hit_pmtid>=300000) { // SPMT
            tool = m_pmt_tool[kPMT_CD_Spmt];
        } else {
            LogWarn << "Unknown pmtid: " << hit_pmtid << std::endl;
            continue;
        }
  
        Pulse pulse = tool->generate(hit_pmtid, hitTimeStamp, hit_amplitude);
        pulse.m_entry=entry;
        pulse.npe=hit->getNPE();
        pulse.hitTime=hit->getHitTime();
        pulse.TTS = hit_timeSpread;
        pulse.timeoffset = m_timeOffset;
        pulse.timeWindow=hit->getTimeWindow();    
        pulse.evtTimeStamp=curTimeStamp;  
        pulse.corrpulseHitTime = corrhitTimeStamp;
        pulse.tof = m_tof;
        if (hit_pmtid<18000) { // LPMT
            LPMT_pulse_buffer.push_back(pulse);
        } else if (hit_pmtid>=300000) { // SPMT
            SPMT_pulse_buffer.push_back(pulse);
        } else {
            LogWarn << "Unknown pmtid: " << hit_pmtid << std::endl;
            continue;
        }
   
        //add afterPulse   
        if(m_enableAfterPulse&&hit_pmtid<=17739){// only LPMT for now
            if(gRandom->Rndm() < m_afterPulseProb){
                double current_rand = gRandom -> Rndm();
                double delta_afterPulseTime = ConvertPdfRand01(current_rand,m_afterPulsePdf,m_afterPulseEdges);
                if(delta_afterPulseTime>10000){
                    delta_afterPulseTime = 10000;
                }
                hitTimeStamp.Add(delta_afterPulseTime * 1e-9);

                if(m_elecsvc->get_isHamamatsu(hit_pmtid)==1){
                    hit_amplitude =  PulseAmp(hit_amplitude,m_gain,m_sigmaGain);
                }else{
                    //TGraph* g1 = m_elecsvc->get_graph();
                    double exp_frac = gNNVTAmp->Eval(0);
                    hit_amplitude =  PulseAmpMCP(hit_amplitude,m_gain,m_sigmaGain,exp_frac);
                }


                Pulse after_pulse = tool->generate(hit_pmtid, hitTimeStamp, hit_amplitude);
                after_pulse.npe=1;
                after_pulse.hitTime = hit->getHitTime() + delta_afterPulseTime;
                after_pulse.timeoffset = m_timeOffset;
                after_pulse.type = kAfterPulse;
                LPMT_pulse_buffer.push_back(after_pulse);
            }
        }
        ++pulses_sum;
    }
    return pulses_sum;

}

int
PMTSimAlg::generate_cd_dark_pulses(TimeStamp starttime,TimeStamp endtime)
{
    std::deque<Pulse>& LPMT_pulse_buffer = m_elecsvc->get_CDLPMT_pulse_buffer();
    std::deque<Pulse>& SPMT_pulse_buffer = m_elecsvc->get_CDSPMT_pulse_buffer();

    int pulses_sum = 0;


    starttime.Subtract(200*1e-9);
 
    TimeStamp delta_Time = endtime - starttime ; 

    double deltaSimTime = delta_Time.GetSeconds();


    int dark_num_total = 0;
    for (int pmtid = 0; pmtid < m_total_Lpmt; ++pmtid) {
        int Ndark = 0;
        if(m_enableAssignDarkPulse) { Ndark = PoissonRand(m_darkRate *deltaSimTime *m_darkRateScale);}                    
        else { Ndark = PoissonRand(( m_elecsvc->get_darkRate(pmtid))*deltaSimTime*m_darkRateScale);}

        for (int dummy = 0; dummy < Ndark; ++dummy) {
            double amplitude = 0.0;
            double darkPulseTime = gRandom->Rndm() * (deltaSimTime*1e9);
            TimeStamp pulseHitTime = starttime;
            pulseHitTime.Add(darkPulseTime*1e-9) ;
            double m_gain = m_elecsvc->get_gain(pmtid);
            double m_sigmaGain = m_elecsvc->get_sigmaGain(pmtid);

            if(m_elecsvc->get_isHamamatsu(pmtid)==1){
                amplitude =  PulseAmp(1.0,m_gain,m_sigmaGain);
            }else{
                //TGraph* g1 = m_elecsvc->get_graph();
                double exp_frac = gNNVTAmp->Eval(0);
                amplitude =  PulseAmpMCP(1.0,m_gain,m_sigmaGain,exp_frac);
            }
 
            IPMTTool* tool = 0;
            tool = m_pmt_tool[kPMT_CD_Lpmt];
            Pulse dark_pulse = tool->generate(pmtid, pulseHitTime, amplitude);
            dark_pulse.npe=1;
            dark_pulse.hitTime = (pulseHitTime-starttime).GetSeconds()*1e9; 
            dark_pulse.corrpulseHitTime = pulseHitTime;
            dark_pulse.type = kDarkPulse;
            dark_pulse.timeoffset = 0;
            dark_pulse.evtTimeStamp = starttime;
            LPMT_pulse_buffer.push_back(dark_pulse);
        }   
        
        dark_num_total += Ndark;
    }

    for (int i = 0; i < m_total_Spmt; ++i) {
        int pmtid = 300000+i;
        int Ndark = PoissonRand(( m_elecsvc->get_darkRate(i+m_total_Lpmt))*deltaSimTime*m_darkRateScale_SPMT); //info need to be in PmtData.root
        for(int dummy = 0; dummy < Ndark; ++dummy){
            double amplitude = 1.0; //1 p.e.
            double darkPulseTime = gRandom->Rndm() * (deltaSimTime*1e9);
            TimeStamp pulseHitTime = starttime;
            pulseHitTime.Add(darkPulseTime*1e-9);

            IPMTTool* tool = 0;
            tool = m_pmt_tool[kPMT_CD_Spmt];
            Pulse dark_pulse = tool->generate(pmtid, pulseHitTime, amplitude);
            dark_pulse.npe=1;
            dark_pulse.hitTime = (pulseHitTime-starttime).GetSeconds()*1e9;
            dark_pulse.corrpulseHitTime = pulseHitTime;
            dark_pulse.type = kDarkPulse;
            dark_pulse.timeoffset = 0;
            dark_pulse.evtTimeStamp = starttime;
            SPMT_pulse_buffer.push_back(dark_pulse);
        }
        dark_num_total += Ndark;
    }
    return dark_num_total;
}
int PMTSimAlg::merge_cd_LPMT_pulses(){
    int pulse_size = 0;
    std::deque<Pulse>& buffer = m_elecsvc->get_CDLPMT_pulse_buffer();
    std::deque<Pulse> buffer_after_merge;
    buffer_after_merge.clear();
    pulse_size = buffer.size();
    LogInfo<<"pulse size before merge:"<<pulse_size<<endl;
    std::deque<Pulse>::iterator it;
    for ( it=buffer.begin(); it != buffer.end(); ++it) {
        if(buffer[1] != buffer[0]){
            buffer_after_merge.push_back(buffer[0]);
        }else{
            buffer[1] += buffer[0];
        }
        buffer.pop_front();
    }
    m_elecsvc->set_CDLPMT_pulse_buffer(buffer_after_merge);
    LogInfo<<"pulse size after merge:"<<buffer_after_merge.size()<<endl;
    pulse_size = buffer_after_merge.size();
    return pulse_size;
}
////////////////////WP pulse///////////////////////////
int PMTSimAlg::generate_wp_pulses(JM::SimEvent* evt, TimeStamp curTimeStamp)
{
    std::deque<Pulse>& WPLPMT_pulse_buffer = m_elecsvc->get_WPLPMT_pulse_buffer();
    int pulses_sum = 0; // generated pulses

    const std::vector<JM::SimPMTHit*>& wp_hits = evt->getWPHitsVec();

    EventKeeper& keeper = EventKeeper::Instance();
    const EventKeeper::Entry& entry = keeper.current_entry();

    for (std::vector<JM::SimPMTHit*>::const_iterator it_hit = wp_hits.begin();
         it_hit != wp_hits.end(); ++it_hit) {

        JM::SimPMTHit* hit=*it_hit;
        int hit_pmtid = hit->getPMTID();
        TimeStamp hitTimeStamp = curTimeStamp;
        double hit_amplitude = 0;

        double m_efficiency = 0;
        double m_afterPulseProb = 0;
        double m_timeSpread = 0;
        double m_timeOffset = 0;
        double m_gain = 0;
        double m_sigmaGain =0;

        if(hit_pmtid>29999&&hit_pmtid<=32400){ // for WP Lpmts
            //LPMT data service
            if(m_enableAssignGain){
                m_gain = m_Gain;
            }else{
                m_gain = m_elecsvc->get_gain(hit_pmtid);
            }
            if(m_enableAssignSigmaGain){
                m_sigmaGain = m_SigmaGain;
            }else{
                m_sigmaGain = m_elecsvc->get_sigmaGain(hit_pmtid);
            }

            m_efficiency = m_elecsvc->get_efficiency(hit_pmtid);
            m_afterPulseProb = m_elecsvc->get_afterPulseProb(hit_pmtid);
            m_timeSpread = m_elecsvc->get_timeSpread(hit_pmtid);
            m_timeOffset = m_elecsvc->get_timeOffset(hit_pmtid);

            // here, the hit time is not the absolute timestamp
            double hit_hitTime = hit->getHitTime() + gRandom->Gaus(0,1) * (m_timeSpread/2.355) + m_timeOffset;
            hitTimeStamp.Add(hit_hitTime*1e-9); // ns -> s  TTS + timeoffset


            double hit_weight = hit->getNPE();
            // add QE affect.
            if(m_enableEfficiency){
                if(gRandom->Rndm() > m_efficiency){
                    continue;    //ignore SimHit due to efficiency
                }
            }

            // single PE resolution
            hit_amplitude =  PulseAmp(hit_weight,m_gain,m_sigmaGain);
        }

        IPMTTool* tool = 0;
        if (hit_pmtid>29999&&hit_pmtid<=32400) { //WP LPMT
            tool = m_pmt_tool[kPMT_WP];
        } else {
            LogWarn << "Unknown pmtid: " << hit_pmtid << std::endl;
            continue;
        }
  
        Pulse pulse = tool->generate(hit_pmtid, hitTimeStamp, hit_amplitude);
        pulse.m_entry=entry;
        // LogWarn << "pulse entry: " << pulse.m_entry.entry << std::endl;
        pulse.npe=hit->getNPE();
        pulse.hitTime=hit->getHitTime();
        pulse.timeWindow=hit->getTimeWindow();    
        pulse.evtTimeStamp=curTimeStamp;
        pulse.npe=hit->getNPE();
        if (hit_pmtid>29999&&hit_pmtid<=32400) { //WP LPMT
            WPLPMT_pulse_buffer.push_back(pulse);
        } else {
            LogWarn << "Unknown pmtid: " << hit_pmtid << std::endl;
            continue;
        }
   
        //add afterPulse   
        if(m_enableAfterPulse&&hit_pmtid>29999&&hit_pmtid<=32400){// WP LPMT for now
            if(gRandom->Rndm() < m_afterPulseProb){
                double current_rand = gRandom -> Rndm();
                double delta_afterPulseTime = ConvertPdfRand01(current_rand,m_afterPulsePdf,m_afterPulseEdges);
                if(delta_afterPulseTime>10000){
                    delta_afterPulseTime = 10000;
                }
                hitTimeStamp.Add(delta_afterPulseTime * 1e-9);
                hit_amplitude = PulseAmp(hit_amplitude,m_gain,m_sigmaGain) ;
                Pulse after_pulse = tool->generate(hit_pmtid, hitTimeStamp, hit_amplitude);
                WPLPMT_pulse_buffer.push_back(after_pulse);
            }
        }
        ++pulses_sum;
    }
    return pulses_sum;

}

int
PMTSimAlg::generate_wp_dark_pulses(TimeStamp starttime,TimeStamp endtime)
{
    std::deque<Pulse>& WPLPMT_pulse_buffer = m_elecsvc->get_WPLPMT_pulse_buffer();

    int pulses_sum = 0;

    starttime.Subtract(200*1e-9);
    endtime.Add(800*1e-9);


    //    LogInfo << "start time : . " <<  starttime <<std::endl;
    //    LogInfo << "end time : . " <<  endtime <<std::endl;

 
    TimeStamp delta_Time = endtime - starttime ; 

    double deltaSimTime = delta_Time.GetSeconds();


    int dark_num_total = 0;
    for (int pmtid = 30000; pmtid <30000+ m_total_WPLpmt; ++pmtid) {    
        int Ndark = 0;
        if(m_enableAssignDarkPulse) { Ndark = PoissonRand(m_darkRate_WP *deltaSimTime* m_darkRateScale_WP);}                    
        else { Ndark = PoissonRand(( m_elecsvc->get_darkRate(pmtid))*deltaSimTime* m_darkRateScale_WP);}

        for (int dummy = 0; dummy < Ndark; ++dummy) {
            double amplitude = 0.0;
            double darkPulseTime = gRandom->Rndm() * (deltaSimTime*1e9);
            TimeStamp pulseHitTime = starttime;
            pulseHitTime.Add(darkPulseTime*1e-9) ;
            double m_gain = m_elecsvc->get_gain(pmtid);
            double m_sigmaGain = m_elecsvc->get_sigmaGain(pmtid);
            amplitude = PulseAmp(1.0, m_gain, m_sigmaGain);
 
            IPMTTool* tool = 0;
            tool = m_pmt_tool[kPMT_WP];
            Pulse dark_pulse = tool->generate(pmtid, pulseHitTime, amplitude);
            WPLPMT_pulse_buffer.push_back(dark_pulse);
        }   
        
        dark_num_total += Ndark;
    }

    return dark_num_total;
}
int PMTSimAlg::merge_wp_LPMT_pulses(){
    int pulse_size = 0;
    std::deque<Pulse>& buffer = m_elecsvc->get_WPLPMT_pulse_buffer();
    std::deque<Pulse> buffer_after_merge;
    buffer_after_merge.clear();
    pulse_size = buffer.size();
    LogInfo<<"pulse size before merge:"<<pulse_size<<endl;
    std::deque<Pulse>::iterator it;
    for ( it=buffer.begin(); it != buffer.end(); ++it) {
        if(buffer[1] != buffer[0]){
            buffer_after_merge.push_back(buffer[0]);
        }else{
            buffer[1] += buffer[0];
        }
        buffer.pop_front();
    }
    m_elecsvc->set_WPLPMT_pulse_buffer(buffer_after_merge);
    LogInfo<<"wp pulse size after merge:"<<buffer_after_merge.size()<<endl;
    pulse_size = buffer_after_merge.size();
    return pulse_size;
}

//////////////////////////////////////

int 
PMTSimAlg::generate_tt_pulses(JM::SimEvent* evt,TimeStamp t)
{
    return 0;
}

std::shared_ptr<JM::EvtNavigator> PMTSimAlg::get_one_event()
{

    // choose one type of data

    // in this prototype, we just load hits from buffer
    // without mixing any other events.

    //    string sample = "default";
    string sample;

    double ranNum = gRandom->Uniform(m_totalRate);
    double sumRate = 0;

    for(map<string,double>::iterator it=m_rateMap.begin(); it!=m_rateMap.end(); it++){
        sumRate += it->second;
        if(ranNum<sumRate){
            sample = it->first;
            break;
        }
    }
    m_current_event.tag = sample;

    if( m_firstMap[sample] ){
        int entries = m_incidentMap[sample]->getEntries();
        //LogInfo<<"initial selected sample: "<<sample <<", sample entries: " <<entries<<endl;

        int num = 0;

        if (m_startIndexMap.count(sample)) {
            num = m_startIndexMap[sample];
        }
        // else {
        //     num = gRandom->Integer(entries-1);
        // }

        //LogInfo << "sample: " << sample <<" first evt starts index: "<<num<<endl;
        m_incidentMap[sample]->reset(num);//we will read data from this entry
        m_firstMap[sample] = false; //for each sample just need set one time
    }


    //LogInfo << " m_incidentMap[" << sample << "]->fire() begin. " << std::endl;

    // In order to handle the missing event, we need to catch the Sniper Incident.
    try {
        m_incidentMap[sample]->fire(*getRoot());
    } catch (StopRunThisEvent& e) {
        return 0;
    } catch (StopRunProcess& e) {
        return 0;
    }
    //LogInfo << " m_incidentMap[" << sample << "]->fire() end. " << std::endl;

    // get current entry
    // m_current_event.entry = m_incidentMap[sample]->getEntry();
    // LogWarn << "current entry: " << m_current_event.entry << std::endl;
    Task* task = m_taskObjMap[sample];
    if (getRoot()->Snoopy().state()==Sniper::RunState::Stopped
        || task->Snoopy().state()==Sniper::RunState::Stopped) {
        LogWarn << "task stopped." << std::endl;
        LogWarn << "top task " << " state: " << static_cast<int>(getRoot()->Snoopy().state()) << std::endl;
        LogWarn << "task " << sample << " state: " << static_cast<int>(task->Snoopy().state()) << std::endl;
        return 0;
    }
    m_current_event.filename = m_incidentMap[sample]->getFileName();

    std::string path = "/Event"; 

    LogDebug << "SniperDataPtr path: " << path << std::endl;
    SniperDataPtr<JM::NavBuffer> navBuf(task,path); 
    if (navBuf.invalid()) {
        LogError << "Can't locate data: " << path << std::endl;
        return 0;
    }
    if (!navBuf->size()) {
        LogWarn << "Can't load more events in " << sample << std::endl;
        return 0;
    }

    if (!navBuf->curEvt()) {
        LogWarn << "Can't load more events in " << sample << std::endl;
        return 0;
    }

    m_current_event.evtnav = *navBuf->current();
    JM::EvtNavigator* Nav = navBuf->curEvt();

    // FIXME
    // To get the real entry number, try to get from SmartRef.
    const std::vector<JM::SmartRef*> refs = Nav->getRef();
    const std::vector<std::string>& paths = Nav->getPath();
    if (refs.size() == 0) {
        LogWarn << "Can't find any SmartRef." << std::endl;
    } else if (refs.size() == 1) {
        const std::string& p = paths[0];
        if (p != "/Event/Sim") {
            LogWarn << "The ref path is not /Event/Sim" << std::endl;
        }

        const JM::SmartRef* ref = refs[0];
        const Long64_t& ref_entry = ref->entry();
        m_current_event.entry = ref_entry;

    } else {
        LogWarn << "There are several SmartRefs in the EvtNavigator."
                << "We will try to find the /Event/Sim."
                << std::endl;

        for (int i = 0; i < refs.size(); ++i) {
            const std::string& p = paths[i];
            if (p != "/Event/Sim") {
                continue;
            }

            const JM::SmartRef* ref = refs[i];
            const Long64_t& ref_entry = ref->entry();

            m_current_event.entry = ref_entry;
        }
    }

    m_current_event.header = dynamic_cast<JM::SimHeader*>(Nav->getHeader("/Event/Sim"));
    //    EventKeeper& keeper = EventKeeper::Instance();
    //    keeper.set_current_entry(m_current_event);
    vec_current_event.push_back(m_current_event); 
    //LogInfo << "CURRENT EVENT: [" << m_current_event.tag << "] "
    //        << m_current_event.filename
    //        << " "
    //        << m_current_event.entry
    //        << std::endl;


    //return Nav;
    return *navBuf->current();
}


std::shared_ptr<JM::EvtNavigator> PMTSimAlg::get_one_SN_event()
{

    string sample = m_SN_sampleName;

    if(sample == "") { LogError << " No SN input file!! " << std::endl; return 0; }

    m_current_event.tag = sample;

    if( m_firstMap[sample] ){
        int entries = m_incidentMap[sample]->getEntries();
        //LogInfo<<"initial selected sample: "<<sample <<", sample entries: " <<entries<<endl;

        int num = 0;

        if (m_startIndexMap.count(sample)) {
            num = m_startIndexMap[sample];
        }
        // else {
        //     num = gRandom->Integer(entries-1);
        // }

        //LogInfo << "sample: " << sample <<" first evt starts index: "<<num<<endl;
        m_incidentMap[sample]->reset(num);//we will read data from this entry
        m_firstMap[sample] = false; //for each sample just need set one time
    }


    //LogInfo << " m_incidentMap[" << sample << "]->fire() begin. " << std::endl;

    // In order to handle the missing event, we need to catch the Sniper Incident.
    try {
        m_incidentMap[sample]->fire(*getRoot());
    } catch (StopRunThisEvent& e) {
        return 0;
    } catch (StopRunProcess& e) {
        return 0;
    }
    //LogInfo << " m_incidentMap[" << sample << "]->fire() end. " << std::endl;

    // get current entry
    // m_current_event.entry = m_incidentMap[sample]->getEntry();
    // LogWarn << "current entry: " << m_current_event.entry << std::endl;
    Task* task = m_taskObjMap[sample];
    if (getRoot()->Snoopy().state()==Sniper::RunState::Stopped
        || task->Snoopy().state()==Sniper::RunState::Stopped) {
        LogWarn << "task stopped." << std::endl;
        LogWarn << "top task " << " state: " << static_cast<int>(getRoot()->Snoopy().state()) << std::endl;
        LogWarn << "task " << sample << " state: " << static_cast<int>(task->Snoopy().state()) << std::endl;
        return 0;
    }
    m_current_event.filename = m_incidentMap[sample]->getFileName();

    std::string path = "/Event"; 

    LogDebug << "SniperDataPtr path: " << path << std::endl;
    SniperDataPtr<JM::NavBuffer> navBuf(task,path); 
    if (navBuf.invalid()) {
        LogError << "Can't locate data: " << path << std::endl;
        return 0;
    }
    if (!navBuf->size()) {
        LogWarn << "Can't load more events in " << sample << std::endl;
        return 0;
    }

    if (!navBuf->curEvt()) {
        LogWarn << "Can't load more events in " << sample << std::endl;
        return 0;
    }

    m_current_event.evtnav = *navBuf->current();
    JM::EvtNavigator* Nav = navBuf->curEvt();

    // FIXME
    // To get the real entry number, try to get from SmartRef.
    const std::vector<JM::SmartRef*> refs = Nav->getRef();
    const std::vector<std::string>& paths = Nav->getPath();
    if (refs.size() == 0) {
        LogWarn << "Can't find any SmartRef." << std::endl;
    } else if (refs.size() == 1) {
        const std::string& p = paths[0];
        if (p != "/Event/Sim") {
            LogWarn << "The ref path is not /Event/Sim" << std::endl;
        }

        const JM::SmartRef* ref = refs[0];
        const Long64_t& ref_entry = ref->entry();
        m_current_event.entry = ref_entry;

    } else {
        LogWarn << "There are several SmartRefs in the EvtNavigator."
                << "We will try to find the /Event/Sim."
                << std::endl;

        for (int i = 0; i < refs.size(); ++i) {
            const std::string& p = paths[i];
            if (p != "/Event/Sim") {
                continue;
            }

            const JM::SmartRef* ref = refs[i];
            const Long64_t& ref_entry = ref->entry();

            m_current_event.entry = ref_entry;
        }
    }

    m_current_event.header = dynamic_cast<JM::SimHeader*>(Nav->getHeader("/Event/Sim"));
    vec_current_event.push_back(m_current_event); 


    //return Nav;
    return *navBuf->current();
}

//cd Lpmt Algorithms
//SPE response
double PMTSimAlg::PulseAmp(double weight,double gain, double sigmaGain){
    double  m_speExpDecay = 1.1;
    double  amp;
    double  m_speCutoff = 0.15; 
    double randW = gRandom->Rndm();
    double m_expWeight = 0.01;
    if (randW > m_expWeight || weight >1.1){
        amp = gRandom->Gaus(0,1) * sigmaGain * TMath::Sqrt(weight) + gain * weight;     
    }   
    else {
        amp = (gRandom->Exp(m_speExpDecay) + m_speCutoff) * gain * weight;
    }   
    if(amp<0) amp = 0;
    
    return amp;
}

double PMTSimAlg::PulseAmpMCP(double weight,double gain, double sigmaGain,double exp_frac){

    double amp;
    double m_speCutoff = 0.1;
    double randW = gRandom->Rndm();
    //enable large signal.
    if (!m_enableMCPLargeSignal ){ randW = 1.01; }
    if(randW > exp_frac){
        amp = gRandom->Gaus(0,1) * sigmaGain * TMath::Sqrt(weight) + gain * weight;
    }
    else {
        amp = (gRandom->Exp(2.2) + m_speCutoff) * gain * weight;
    }
    if(amp<0) amp = 0;
    return amp;
}

double PMTSimAlg::ConvertPdfRand01(double rand,vector<double> pdf, vector<double> edges){
    // Defined PDF returns random number in [0,1] distributed according to user-defined histogram.
    // It assumes even bin sizes, so accomodate uneven bin sizes for generality. 
    int current_bin = 0;
    int Nbins = pdf.size();

    for(int bin=0; bin<Nbins; bin++) {
        if(rand >= pdf[bin] && rand < pdf[bin+1]) {
            current_bin = bin;
            break;
        }
        else
            current_bin = Nbins-1;
    }


    return edges[current_bin] + (rand-pdf[current_bin])*(edges[current_bin+1]-edges[current_bin])
        /(pdf[current_bin+1]-pdf[current_bin]);

}    

int PMTSimAlg::PoissonRand(double mean) {

    int n;
    if (mean <= 0) return 0;

    double expmean = exp(-mean);
    double pir = 1;
    n = -1;
    while(1) {
        n++;
        pir *= gRandom->Rndm();
        if (pir <= expmean) break;
    }
    return n;
}


