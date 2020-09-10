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

    dyn_pde_old = np.array(dyn_pde_old) * 0.872
    hmcp_pde_old = np.array(hmcp_pde_old) * 0.957
    lmcp_pde_old = np.array(lmcp_pde_old) * 0.916
    all_pde_old = list(dyn_pde_old) + list(hmcp_pde_old) + list(lmcp_pde_old); 

    all_pde_new, dyn_pde_new, hmcp_pde_new, lmcp_pde_new = [], [], [], []
    for index, ham_flag, pde_flag, pde2 in zip(pmtid, isHam, isHiQE, pde_new):
        if ham_flag :
            dyn_pde_new.append(pde2*0.98)
            #all_pde_new.append(pde2*0.98)
            all_pde_new.append(pde2*0.98)
        elif not ham_flag and pde_flag:
            hmcp_pde_new.append(pde2*0.96)
            #all_pde_new.append(pde2*0.96)
            all_pde_new.append(pde2*0.96)
        elif not ham_flag and not pde_flag:
            lmcp_pde_new.append(pde2*0.96)
            all_pde_new.append(pde2*0.96)
            #all_pde_new.append(pde2*0.96)
        else:
            print("Unknow PMT Type ??")


    dyn_qe_old, dyn_qe_new = [], []
    hmcp_qe_old, hmcp_qe_new = [], []
    lmcp_qe_old, lmcp_qe_new = [], []
    all_qe_old, all_qe_new = [], []
    dyn_qe_old = np.array(dyn_pde_old) / 0.872
    dyn_qe_new = np.array(dyn_pde_new) / 0.872
    hmcp_qe_old = np.array(hmcp_pde_old) / 0.916
    hmcp_qe_new = np.array(hmcp_pde_new) / 0.916
    lmcp_qe_old = np.array(lmcp_pde_old) / 0.957
    lmcp_qe_new = np.array(lmcp_pde_new) / 0.957
    
    all_qe_old = list(dyn_qe_old) + list(hmcp_qe_old) + list(lmcp_qe_old)
    all_qe_new = list(dyn_qe_new) + list(hmcp_qe_new) + list(lmcp_qe_new)

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
    """

    plt.figure(0)
    plt.hist(dyn_qe_old, bins=100, range=(20, 40), alpha=0.4, color="salmon", label="J19 HamPmt: %.2f"%np.array(dyn_qe_old).mean())
    plt.hist(dyn_qe_new, bins=100, range=(20, 40), alpha=0.4, color="cornflowerblue", label="J20 HamPmt: %.2f"%np.array(dyn_qe_new).mean())
    print("dyn old pmt number: %d" %len(dyn_qe_old))
    print("dyn new pmt number: %d" %len(dyn_qe_new))
    plt.legend()
    plt.xlabel("QE")
    plt.title("QE")
    #plt.show()
    plt.savefig("DynPmt_QE.pdf")

    plt.figure(1)
    plt.hist(hmcp_qe_old, bins=100, range=(20, 40), alpha=0.4, color="salmon", label="J19 HiQE MCPPmt: %.2f"%np.array(hmcp_qe_old).mean())
    plt.hist(hmcp_qe_new, bins=100, range=(20, 40), alpha=0.4, color="cornflowerblue", label="J20 HiQE MCPPmt: %.2f"%np.array(hmcp_qe_new).mean())
    print("mcp old pmt number: %d" %len(hmcp_qe_old))
    print("mcp new pmt number: %d" %len(hmcp_qe_new))
    plt.legend()
    plt.xlabel("QE")
    plt.title("QE")
    #plt.show()
    plt.savefig("HMCPPmt_QE.pdf")


    plt.figure(2)
    plt.hist(lmcp_qe_old, bins=100, range=(20, 40), alpha=0.4, color="salmon", label="J19 Normal MCPPmt: %.2f"%np.array(lmcp_qe_old).mean())
    plt.hist(lmcp_qe_new, bins=100, range=(20, 40), alpha=0.4, color="cornflowerblue", label="J20 Normal MCPPmt: %.2f"%np.array(lmcp_qe_new).mean())
    print("all old pmt number: %d" %len(lmcp_qe_old))
    print("all new pmt number: %d" %len(lmcp_qe_new))
    plt.legend()
    plt.xlabel("QE")
    plt.title("QE")
    #plt.show()
    plt.savefig("LMCPPmt_QE.pdf")

    plt.figure(3)
    plt.hist(all_qe_old, bins=100, range=(20, 40), alpha=0.4, color="salmon", label="J19 AllPmt: %.2f"%np.array(all_qe_old).mean())
    plt.hist(all_qe_new, bins=100, range=(20, 40), alpha=0.4, color="cornflowerblue", label="J20 AllPmt: %.2f"%np.array(all_qe_new).mean())
    print("all old pmt number: %d" %len(all_qe_old))
    print("all new pmt number: %d" %len(all_qe_new))
    plt.legend()
    plt.xlabel("QE")
    plt.title("QE")
    #plt.show()
    plt.savefig("AllPmt_QE.pdf")

    

