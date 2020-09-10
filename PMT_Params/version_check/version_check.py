#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Tue Sep  8 12:13:12 2020
# File Name: version_check.py
"""

import matplotlib.pyplot as plt
import numpy as np
from ROOT import TFile, TTree
import uproot as up

def read_J20():
    infile = up.open("../Other_Params/PmtData_Lpmt0908.root")
    tree = infile["PmtData_Lpmt"]
    pmtid = tree["pmtID"].array()
    isHam = tree["MCP_Hama"].array()
    isHiQE = tree["HiQE_MCP"].array()
    pde = tree["PDE"].array()
    dcr = tree["DCR"].array()
    gain = tree["Gain"].array() / 1e7
    sigma = tree["Resolution"].array()
    tts = tree["TTS_SS"].array()
    return pmtid, isHam, isHiQE, pde, dcr, gain, sigma, tts


def read_oldPde():
    dyn_pde_old, hmcp_pde_old, lmcp_pde_old = [], [], []
    with open("/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v1r0-Pre2/data/Simulation/ElecSim/pmtdata.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            if data[1]=="Hamamatsu":
                dyn_pde_old.append(float(data[2])*100)
            elif data[1] == "HighQENNVT":
                hmcp_pde_old.append(float(data[2])*100)
            elif data[1] == "NNVT":
                lmcp_pde_old.append(float(data[2])*100)
    return dyn_pde_old, hmcp_pde_old, lmcp_pde_old
        

def read_oldroot():
    infile = up.open("/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v1r0-Pre2/data/Simulation/ElecSim/PmtData.root")
    tree = infile["PmtData"]
    dcr  = tree["darkRate"].array() / 1000
    gain = tree["gain"].array()
    sigma = tree["sigmaGain"].array()
    tts = tree["timeSpread"].array() / 2.35
    return dcr, gain, sigma, tts

    


#plt.style.use("seaborn-paper")


if __name__ == "__main__":
    pmtid, isHam, isHiQE, pde_new, dcr_new, gain_new, rsl_new, tts_new = read_J20()
    dyn_pde_old, hmcp_pde_old, lmcp_pde_old = read_oldPde()
    dcr_old, gain_old, rsl_old, tts_old = read_oldroot()

    """
    all_pde_old = dyn_pde_old + hmcp_pde_old + lmcp_pde_old; 

    all_pde_new, dyn_pde_new, hmcp_pde_new, lmcp_pde_new = [], [], [], []
    for index, ham_flag, qe_flag, pde2 in zip(pmtid, isHam, isHiQE, pde_new):
        if ham_flag :
            dyn_pde_new.append(pde2/0.872*0.98)
            #all_pde_new.append(pde2*0.98)
            all_pde_new.append(pde2/0.872*0.98)
            print("%d  %.3f  %.3f" %(index, pde2*0.98, all_pde_new[-1]) )
        elif not ham_flag and qe_flag:
            hmcp_pde_new.append(pde2/0.916*0.96)
            #all_pde_new.append(pde2*0.96)
            all_pde_new.append(pde2/0.916*0.96)
            print("%d  %.3f  %.3f" %(index, pde2*0.96, all_pde_new[-1]) )
        elif not ham_flag and not qe_flag:
            lmcp_pde_new.append(pde2/0.916*0.96)
            all_pde_new.append(pde2/0.916*0.96)
            #all_pde_new.append(pde2*0.96)
            print("%d  %.3f  %.3f" %(index, pde2*0.96, all_pde_new[-1]) )
        else:
            print("Unknow PMT Type ??")

    #for idx, qe in enumerate(all_pde_new):
    #    print("%d %.2f"%(idx, qe))

    """

    dyn_tts_old, dyn_tts_new = [], []
    mcp_tts_old, mcp_tts_new = [], []
    all_tts_old, all_tts_new = [], []
    for ham_flag, oldtts, newtts in zip(isHam, tts_old, tts_new):
        if ham_flag:
            dyn_tts_old.append(oldtts)
            all_tts_old.append(oldtts)
            dyn_tts_new.append(newtts)
            all_tts_new.append(newtts)
        else:
            mcp_tts_old.append(oldtts)
            all_tts_old.append(oldtts)
            mcp_tts_new.append(newtts)
            all_tts_new.append(newtts)


    plt.figure(0)
    plt.hist(dyn_tts_old, bins=100, range=(0, 4), alpha=0.4, color="salmon", label="J19 AllPmt: %.2f"%np.array(dyn_tts_old).mean())
    plt.hist(dyn_tts_new, bins=100, range=(0, 4), alpha=0.4, color="cornflowerblue", label="J20 AllPmt: %.2f"%np.array(dyn_tts_new).mean())
    print("dyn old pmt number: %d" %len(dyn_tts_old))
    print("dyn new pmt number: %d" %len(dyn_tts_new))
    plt.legend()
    plt.xlabel("TTS/ns")
    plt.title("TTS")
    #plt.show()
    plt.savefig("DynPmt_TTS.pdf")

    plt.figure(1)
    plt.hist(mcp_tts_old, bins=100, range=(0, 10), alpha=0.4, color="salmon", label="J19 AllPmt: %.2f"%np.array(mcp_tts_old).mean())
    plt.hist(mcp_tts_new, bins=100, range=(0, 10), alpha=0.4, color="cornflowerblue", label="J20 AllPmt: %.2f"%np.array(mcp_tts_new).mean())
    print("mcp old pmt number: %d" %len(mcp_tts_old))
    print("mcp new pmt number: %d" %len(mcp_tts_new))
    plt.legend()
    plt.xlabel("TTS/ns")
    plt.title("TTS")
    #plt.show()
    plt.savefig("MCPPmt_TTS.pdf")

    """
    plt.figure(2)
    plt.hist(all_tts_old, bins=100, range=(0.5, 1.5), alpha=0.4, color="salmon", label="J19 AllPmt: %.2f"%np.array(all_tts_old).mean())
    plt.hist(all_tts_new, bins=100, range=(0.5, 1.5), alpha=0.4, color="cornflowerblue", label="J20 AllPmt: %.2f"%np.array(all_tts_new).mean())
    print("all old pmt number: %d" %len(all_tts_old))
    print("all new pmt number: %d" %len(all_tts_new))
    plt.legend()
    plt.xlabel("Gain")
    plt.title("Gain")
    #plt.show()
    plt.savefig("AllPmt_Gain.pdf")
    """


    

