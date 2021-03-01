from ROOT import TH2D, TFile

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

import numpy as np
def tof_calc(sx, sy, sz, px, py, pz):
    rindex = 1.4931
    c = 299.792458
    v = c/rindex
    d = np.sqrt( (sx-px)**2 + (sy-py)**2 + (sz-pz)**2 )
    return d/v


def loadCubic():
    cubicx, cubicy, cubicz = [], [], []
    with open("cubic_center.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            cubicx.append(float(data[0]))
            cubicy.append(float(data[1]))
            cubicz.append(float(data[2]))
    return cubicx, cubicy, cubicz
    

import matplotlib.pyplot as plt
def main():
    pmtx, pmty, pmtz = loadPmtPos()
    cubicx, cubicy, cubicz = loadCubic()

    hist = TH2D("tof_map", "", 17612, 0, 17612, 179, 0, 179)
    #tof_map = [[] for i in range(179)]
    for i in range(179):
        print("Cubic volume No. %d" %i )
        sx, sy, sz = cubicx[i], cubicy[i], cubicz[i]
        for j in range(17612):
            px, py, pz = pmtx[j], pmty[j], pmtz[j]
            tt = tof_calc(sx, sy, sz, px, py, pz)
            #tof_map[i].append(tt)
            hist.SetBinContent(j+1, i+1, tt)

    hist.SaveAs("tof_map.root")
    #print(tof_map[0])
    #plt.plot(tof_map[0], "-")
    #plt.imshow(tof_map, interpolation='nearest', origin='lower', \
    #           extent=[0, 17612, 0, 179])
    #plt.show()



if __name__ == "__main__":
    main()
