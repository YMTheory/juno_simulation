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
    return pmtid, isHam, isHiQE, pde, dcr


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
    dcr = tree["darkRate"].array()

    


#plt.style.use("seaborn-paper")


if __name__ == "__main__":
    pmtid, isHam, isHiQE, pde_new, dcr_new = read_J20()
    dyn_pde_old, hmcp_pde_old, lmcp_pde_old = read_oldPde()

    all_pde_old = dyn_pde_old + hmcp_pde_old + lmcp_pde_old; 

    all_pde_new, dyn_pde_new, hmcp_pde_new, lmcp_pde_new = [], [], [], []
    for index, ham_flag, qe_flag, pde2 in zip(pmtid, isHam, isHiQE, pde_new):
        if ham_flag :
            dyn_pde_new.append(pde2/0.872*0.98)
            all_pde_new.append(pde2/0.872*0.98)
        elif not ham_flag and qe_flag:
            hmcp_pde_new.append(pde2/0.916*0.96)
            all_pde_new.append(pde2/0.916*0.96)
        elif not ham_flag and not qe_flag:
            lmcp_pde_new.append(pde2/0.916*0.96)
            all_pde_new.append(pde2/0.916*0.96)
        else:
            print("Unknow PMT Type ??")


    plt.hist(all_pde_old, bins=100, range=(20, 40), alpha=0.4, color="salmon", label="J19 AllPmt: %.2f"%np.array(all_pde_old).mean())
    plt.hist(all_pde_new, bins=100, range=(20, 40), alpha=0.4, color="cornflowerblue", label="J20 AllPmt: %.2f"%np.array(all_pde_new).mean())
    print("all old pmt number: %d" %len(all_pde_old))
    print("all new pmt number: %d" %len(all_pde_new))
    plt.legend()
    plt.xlabel("QE")
    plt.title("QE@420nm")
    #plt.show()
    plt.savefig("AllPmt_QE.pdf")



