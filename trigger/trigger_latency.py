import uproot as up
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt

def loaddata(filename):
    ff = up.open(filename)
    pmtid = ff["evt"].array("pmtID")
    hittime = ff["evt"].array("hitTime")

    return pmtid, hittime


def loadPmtPos():
    pmtx, pmty, pmtz = [], [], []
    with open("PMTPos_Acrylic_with_chimney.csv") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            pmtx.append(float(data[1]))
            pmty.append(float(data[2]))
            pmtz.append(float(data[3]))
    return pmtx, pmty, pmtz

pmtx, pmty, pmtz = loadPmtPos()

def tof_calc(sx, sy, sz, px, py, pz):
    rindex = 1.4931
    c = 299.792458
    v = c/rindex
    d = np.sqrt( (sx-px)**2 + (sy-py)**2 + (sz-pz)**2 )
    return d/v



def pulse_buffer(pmtid, hittime, corr):
    if not corr:
        data_array = np.vstack([pmtid, hittime]).T

    else:  # do tof correction
        for i, j in zip(pmtid, hittime):
            if i<20000:
                j -= tof_calc(0, 0, 0, pmtx[int(i)], pmty[int(i)], pmtz[int(i)])
        data_array = np.vstack([pmtid, hittime]).T

    df = pd.DataFrame(data_array, index=None, columns=['pmtid', 'hittime'])
    df.sort_values("hittime", inplace=True)
    return df




trigger_window = 48 #ns
pmt_threshold = 100
trigger_slip = 16 #ns

readout_window = 1000 #ns



def trigger_search(df, beginTime, endTime):
    FiredPMT = [ 0 for i in range(17612) ]
    found = False

    arr = df.to_numpy().T
    for i, j in zip(arr[0], arr[1]):
        if i<20000 and beginTime<j<endTime:
            FiredPMT[int(i)] = 1
        #else:
        #    break

    totFired = 0
    for ii in FiredPMT:
        if ii:
            totFired += 1

    if totFired >= pmt_threshold:
        #print("Find a trigger at : %d ns with %d PMT fired" %(beginTime, totFired) )
        found = True
    else:
        # move trigger window:
        pass

    return found, beginTime
        
        

def slip_window(maxTime, df):
    trigger_buffer = []
    beginTime, endTime = 0, 0
    while beginTime < maxTime:
        endTime += trigger_window
        #print("trigger window range: [%d, %d]" %(beginTime, endTime) )
        found, time = trigger_search(df, beginTime, endTime)
        if found:
            trigger_buffer.append(time)
        beginTime = beginTime + trigger_slip
        endTime = beginTime
    
    new_trigger_buffer = []
    if len(trigger_buffer) > 0:
        atrigger = trigger_buffer[0]
        new_trigger_buffer.append(atrigger)
        for i in range(len(trigger_buffer)-1):
            if trigger_buffer[i+1] > atrigger + readout_window: # out readout window
                atrigger = trigger_buffer[i+1]
                new_trigger_buffer.append(atrigger)

    return new_trigger_buffer

def loadpmt():
    ff = up.open("./PmtData_Lpmt.root")
    tts  = ff["PmtData_Lpmt"].array("TTS_SS")
    toff = ff["PmtData_Lpmt"].array("timeOffset")
    dcr  = ff["PmtData_Lpmt"].array("DCR")
    return tts, toff, dcr

tts, toff, dcr = loadpmt()


def pmt_time(pmtid):
    tcorr = toff[pmtid] + np.random.normal(0, tts[pmtid])
    return tcorr


def add_dcr(beginTime, endTime):
    pmtid_dn, hittime_dn = [], []
    delta = endTime - beginTime
    for i in range(17612):
        num = int( np.random.poisson(dcr[i]*1000*delta/1e9) )
        for j in range(num):
            t = np.random.uniform(beginTime, endTime)
            pmtid_dn.append(i)
            hittime_dn.append(t)

    return pmtid_dn, hittime_dn

        
def toy_time(pmtid, hittime):
    for i, j in zip(pmtid, hittime):
        if i<20000:
            j += pmt_time(int(i))

    return pmtid, hittime




def main():
    filename = "/junofs/production/data-production/Pre-Releases/J20v2r0-Pre0/ACU+CLS/Laser/photon_11522/Laser_0_0_0/detsim/user-root/user-detsim-80013.root"
    pmtid, hittime = loaddata(filename)
    
    trigger_container = []

    evtid = 0
    for arr0, arr1 in zip(pmtid, hittime):
        print("EVENTID -> %d" %evtid)

        minTime = 0
        maxTime = arr1.max() + 300

        pmtid_dn, hittime_dn = add_dcr(minTime, maxTime)
        #print("Total DN geneated : %d" %len(hittime_dn))

        arr2 = np.hstack([arr0, pmtid_dn])
        arr3 = np.hstack([arr1, hittime_dn])

        arr4, arr5 = toy_time(arr2, arr3)

        df = pulse_buffer(arr4, arr5, True)
        trigger_buffer = slip_window(maxTime, df)

        #print(trigger_buffer)
        for tg in trigger_buffer:
            trigger_container.append(tg)

        evtid += 1

    led_time = np.random.uniform(-8, 8, len(trigger_container) )
    plt.hist(trigger_container-led_time, bins=40, alpha=0.7)
    plt.xlabel("latency/ns")
    plt.show()




if __name__ == "__main__":
    main()
