#!/usr/bin/env python
# coding=utf-8

import pandas as pd
from ROOT import TH2D, TH1D, TFile, TRandom, TTree
import ROOT
import h5py
import uproot as up
import pdb
import numpy as np
import matplotlib.pyplot as plt
import random


def read_file():
    filename = "./Pan-Asia-container-barePMT-passed-limin.xlsx"
    df = pd.read_excel(filename, sheet_name="data")
    return df

def pde_dcr_analysis(df, pmttype, isHqe):   # read test data from xlsx from ZhongShan
    """ pmttyp --> MCP 1 , Ham: 0"""
    pde = df["PDE"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ]
    dcr = df["DCR"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ]
    print("Current type pmt: %d" %len(pde) )
    return pde, dcr
    

def read_alldata(df, pmttype, isHqe):
    sn = df["SN"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    gain = df["Gain"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    rsl = df["Resolution"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    tts = df["TTS"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    pde = df["PDE"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    dcr = df["DCR"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    ishqe = df["HiQE_MCP"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    ismcp = df["MCP_Hama"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    amp = df["Amplitude"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    hv = df["HV"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    pvsv = df["PvsV"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    svsn = df["SvsN"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    riset = df["Risetime"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    fallt = df["Falltime"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    fwhm = df["FWHM"].loc[(df["MCP_Hama"]==pmttype) & (df["HiQE_MCP"]==isHqe) ].to_numpy()
    

    return sn, gain, rsl, pde, dcr, tts, amp, hv, pvsv, svsn, riset, fallt, fwhm




def sort_lqmcp(df):
    sn3, gain3, rsl3, pde3, dcr3, tts3, amp3, hv3, pvsv3, svsn3, riset3, fallt3, fwhm3 = read_alldata(df, 1, False)
    # prepare faked data for non-provided pmts
    """ LQMCP only 2032 highest PDE ones used """
    index = [i for i in range(len(pde3))]
    index_sort = np.array( [ i for _,i in sorted(zip(pde3, index), reverse=True)] )
    
    sn3_sorted, gain3_sorted, rsl3_sorted, pde3_sorted, dcr3_sorted = [], [], [], [], []
    tts3_sorted, amp3_sorted, hv3_sorted, pvsv3_sorted, svsn3_sorted, riset3_sorted, fallt3_sorted, fwhm3_sorted = [], [], [], [], [], [], [], []
    for i in index_sort[0:2032]:
        sn3_sorted.append(sn3[i])
        gain3_sorted.append(gain3[i])
        rsl3_sorted.append(rsl3[i])
        dcr3_sorted.append(dcr3[i])
        pde3_sorted.append(pde3[i])
        tts3_sorted.append(tts3[i])
        amp3_sorted.append(amp3[i])
        hv3_sorted.append(hv3[i])
        pvsv3_sorted.append(pvsv3[i])
        svsn3_sorted.append(svsn3[i])
        riset3_sorted.append(riset3[i])
        fallt3_sorted.append(fallt3[i])
        fwhm3_sorted.append(fwhm3[i])

    print("Final selected lqmcp number : %d" % len(sn3_sorted) )

    return sn3_sorted, gain3_sorted, rsl3_sorted, pde3_sorted, dcr3_sorted, tts3_sorted, amp3_sorted, hv3_sorted, pvsv3_sorted, svsn3_sorted, riset3_sorted, fallt3_sorted, fwhm3_sorted


def sample_pmt(df, ismcp, ishq, num):
    sn2, gain2, rsl2, pde2, dcr2, tts2, amp2, hv2, pvsv2, svsn2, riset2, fallt2, fwhm2 = read_alldata(df, ismcp, ishq)
    """ HQMCP require 3456 more added """
    """ DynodeMCP require 5 more added """
    """ Prepare sampling hists """
    hqmcp_gainrsl_hist = TH2D("hqmcp_gainrsl_hist", "", 100, 9000000, 13000000, 30, 0.2, 0.5)
    hqmcp_pdedcr_hist = TH2D("hqmcp_pdedcr_hist", "", 100, 0, 50, 100, 0, 100)
    hqmcp_tts_hist = TH1D("tts", "", 100, 0, 50)
    hqmcp_amp_hist = TH1D("amp", "", 100, 0, 50)
    hqmcp_pvsv_hist = TH1D("pvsv", "", 100, 0, 30)
    hqmcp_svsn_hist = TH1D("svsn", "", 100, 0, 0.2)
    hqmcp_rise_hist = TH1D("rise", "", 100, 0, 50)
    hqmcp_fall_hist = TH1D("fall", "", 100, 0, 50)
    hqmcp_fwhm_hist = TH1D("fwhm", "", 100, 0, 50)
    hqmcp_hv_hist = TH1D("hv", "", 100, 0, 3000)
    for i, j, k, l, a, b, c, d, e, f, g, h in zip(gain2, rsl2, pde2, dcr2, tts2, amp2, hv2, pvsv2, svsn2, riset2, fallt2, fwhm2):
        hqmcp_gainrsl_hist.Fill(i, j)
        hqmcp_pdedcr_hist.Fill(k, l)
        hqmcp_tts_hist.Fill(a)
        hqmcp_amp_hist.Fill(b)
        hqmcp_hv_hist.Fill(c)
        hqmcp_pvsv_hist.Fill(d)
        hqmcp_svsn_hist.Fill(e)
        hqmcp_rise_hist.Fill(f)
        hqmcp_fall_hist.Fill(g)
        hqmcp_fwhm_hist.Fill(h)

    gain2_add, rsl2_add, pde2_add, dcr2_add = [], [], [], []
    tts2_add, amp2_add, hv2_add, pvsv2_add, svsn2_add, riset2_add, fallt2_add, fwhm2_add = [], [], [], [], [], [], [], []
    for idx in range(num):
        pde_sample, dcr_sample = ROOT.Double(0), ROOT.Double(0)
        hqmcp_pdedcr_hist.GetRandom2(pde_sample, dcr_sample)
        gain_sample, rsl_sample = ROOT.Double(0), ROOT.Double(0)
        hqmcp_gainrsl_hist.GetRandom2(gain_sample, rsl_sample)
        gain2_add.append(gain_sample)
        rsl2_add.append(rsl_sample)
        pde2_add.append(pde_sample)
        dcr2_add.append(dcr_sample)
        tts2_add.append(hqmcp_tts_hist.GetRandom())
        amp2_add.append(hqmcp_amp_hist.GetRandom())
        hv2_add.append(hqmcp_hv_hist.GetRandom())
        pvsv2_add.append(hqmcp_pvsv_hist.GetRandom())
        svsn2_add.append(hqmcp_svsn_hist.GetRandom())
        riset2_add.append(hqmcp_rise_hist.GetRandom())
        fallt2_add.append(hqmcp_fall_hist.GetRandom())
        fwhm2_add.append(hqmcp_fwhm_hist.GetRandom())

    gain2_add = np.array(gain2_add) ; gain2_all = np.concatenate((gain2, gain2_add))
    rsl2_add = np.array(rsl2_add) ; rsl2_all = np.concatenate((rsl2, rsl2_add))
    pde2_add = np.array(pde2_add) ; pde2_all = np.concatenate((pde2, pde2_add))
    dcr2_add = np.array(dcr2_add) ; dcr2_all = np.concatenate((dcr2, dcr2_add))
    tts2_add = np.array(tts2_add) ; tts2_all = np.concatenate((tts2, tts2_add))
    amp2_add = np.array(amp2_add) ; amp2_all = np.concatenate((amp2, amp2_add))
    hv2_add = np.array(hv2_add) ; hv2_all = np.concatenate((hv2, hv2_add))
    pvsv2_add = np.array(pvsv2_add) ; pvsv2_all = np.concatenate((pvsv2, pvsv2_add))
    svsn2_add = np.array(svsn2_add) ; svsn2_all = np.concatenate((svsn2, svsn2_add))
    riset2_add = np.array(riset2_add) ; riset2_all = np.concatenate((riset2, riset2_add))
    fallt2_add = np.array(fallt2_add) ; fallt2_all = np.concatenate((fallt2, fallt2_add))
    fwhm2_add = np.array(fwhm2_add) ; fwhm2_all = np.concatenate((fwhm2, fwhm2_add))
            
    return  sn2, gain2_all, rsl2_all, pde2_all, dcr2_all, tts2_all, amp2_all, hv2_all, pvsv2_all, svsn2_all, riset2_all, fallt2_all, fwhm2_all


def faked_sn(pmttype):
    """ faked PMT SN: 0 for Ham, 1 for HqMCP """
    if pmttype == 1:
        hqmcp_sn_faked = ["PA{:d}F".format(i) for i in range(3456) ]
        return hqmcp_sn_faked
    elif pmttype == 0 :
        dyn_sn_faked = ["EA{:d}F" .format(i) for i in range(5) ]
        return dyn_sn_faked



def read_dynID():
    dyn_id = []
    with open("./Hamamatsu_pmtID.txt") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            dyn_id.append(int(line))
    return dyn_id


def gen_lqmcpID(dyn_id):
    num, lqmcp_id = 0, []
    while num<2032:
        sample_id = int( random.uniform(0, 17612) )
        if sample_id in dyn_id or sample_id in lqmcp_id:
            continue;
        lqmcp_id.append(sample_id)
        num += 1

    #lqmcp_id = sorted(lqmcp_id)
    return lqmcp_id
        
def gen_hqmcpID(dyn_id, lqmcp_id):
    num, hqmcp_id = 0, []
    for i in range(17612):
        if i in dyn_id: 
            continue
        elif i in lqmcp_id:
            continue
        hqmcp_id.append(i)
    #hqmcp_id = sorted(hqmcp_id)
    return hqmcp_id


def read_pmtPos():
    df = pd.read_csv("PMTPos_Acrylic_with_chimney.csv", sep=" ")
    return df["x"].to_numpy(), df["y"].to_numpy(), df["z"].to_numpy()



if __name__ == "__main__" :


    """
    # correlation analysis
    pde, dcr = pde_dcr_analysis(read_file(), 1, False)
    hist2d = TH2D("lqmcp2d", "", 100, 0, 50, 100, 0, 100);  
    histpde = TH1D("lqmcp_pde", "", 100, 0, 50)
    
    istdcr = TH1D("lqmcp_dcr", "", 100, 0, 100)
    for p, d in zip(pde ,dcr):
        hist2d.Fill(p, d)
        histpde.Fill(p)
        histdcr.Fill(d)

    out = TFile("pde2dcr_lqmcp.root", "recreate")
    hist2d.Write()
    histpde.Write()
    histdcr.Write()
    out.Close()
    """

    plt.style.use("bmh")
    

    """ read from xlxs --> new test data """
    df = read_file()
    
    # dynode
    sn1_fake = faked_sn(0) 
    dyn_id = read_dynID(); print(len(dyn_id))
    #gain1, rsl1, pde1, dcr1 = sample_pmt(df, 0, False, 5)
    sn1_origin, gain1, rsl1, pde1, dcr1, tts1, amp1, hv1, pvsv1, svsn1, riset1, fallt1, fwhm1 = sample_pmt(df, 0, False, 5)
    sn1 = np.concatenate((sn1_origin, np.array(sn1_fake)) )
    isDyn1 = [1 for i in range(len(dyn_id))]
    isHqe1 = [0 for i in range(len(dyn_id))]


    # lqmcp data
    #sn1, gain1, rsl1, pde1, dcr1 = sort_lqmcp(df)
    sn2, gain2, rsl2, pde2, dcr2, tts2, amp2, hv2, pvsv2, svsn2, riset2, fallt2, fwhm2 = sort_lqmcp(df)
    lqmcp_id = gen_lqmcpID(dyn_id); print(len(lqmcp_id))
    isDyn2 = [0 for i in range(len(lqmcp_id ))]
    isHqe2 = [0 for i in range(len(lqmcp_id)) ]


    # hqmcp
    sn3_fake = faked_sn(1) 
    #gain3, rsl3, pde3, dcr3 = sample_pmt(df, 1, True, 3456)
    sn3_origin, gain3, rsl3, pde3, dcr3, tts3, amp3, hv3, pvsv3, svsn3, riset3, fallt3, fwhm3 = sample_pmt(df, 1, True, 3456)
    sn3 = np.concatenate((sn3_origin, np.array(sn3_fake)) )
    pmt_id = [i for i in range(17612) ]
    hqmcp_id = gen_hqmcpID(dyn_id, lqmcp_id); print(len(hqmcp_id))
    isDyn3 = [0 for i in range(len(hqmcp_id ))]
    isHqe3 = [1 for i in range(len(hqmcp_id)) ]

    tts_ss_dyn = []
    with open("dynode_tts_sampled.csv") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            tts_ss_dyn.append(float(line))
    

    tts_ss_mcp, tts_ss_lqmcp, tts_ss_hqmcp = [], [], []
    with open("mcp_tts_sampled.csv") as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            tts_ss_mcp.append(float(line))

    infile = up.open("PmtData.root") 
    toff = infile["PmtData"].array("timeOffset")
    app = infile["PmtData"].array("afterPulseProb")
    ppp = infile["PmtData"].array('prePulseProb')

    
    idx_final, isDyn_final, isHqe_final, sn_final = [], [], [], []
    gain_final, rsl_final, tts_final, pde_final = [], [], [], []
    dcr_final, hv_final, pvsv_final, svsn_final = [], [], [], []
    amp_final, riset_final, fallt_final, fwhm_final = [], [], [], []
    tts_ss_final, toff_final, app_final, ppp_final = [], [], [], []
    pmtx_final, pmty_final, pmtz_final = read_pmtPos()


    num = 0
    print(np.array(pde2).shape)
    for i in range(17612):
        toff_final.append(toff[i]); ppp_final.append(ppp[i]); app_final.append(app[i])
        if i in dyn_id:
            idx = dyn_id.index(i)
            idx_final.append(dyn_id[idx])
            #print("{:d}, {:d}, {:.2f}" .format(i, idx, pde1[idx]))
            sn_final.append(sn1[idx]); gain_final.append(gain1[idx]); rsl_final.append(rsl1[idx]); dcr_final.append(dcr1[idx])
            pde_final.append(pde1[idx]); isDyn_final.append(isDyn1[idx]); isHqe_final.append(isHqe1[idx]); tts_final.append(tts1[idx])
            amp_final.append(amp1[idx]); hv_final.append(hv1[idx]); pvsv_final.append(pvsv1[idx]); svsn_final.append(svsn1[idx])
            riset_final.append(riset1[idx]); fallt_final.append(fallt1[idx]); fwhm_final.append(fwhm1[idx]); tts_ss_final.append(tts_ss_dyn[idx])
            num += 1
        elif i in lqmcp_id:
            idx = lqmcp_id.index(i)
            idx_final.append(lqmcp_id[idx])
            sn_final.append(sn2[idx]); gain_final.append(gain2[idx]); rsl_final.append(rsl2[idx]); dcr_final.append(dcr2[idx])
            pde_final.append(pde2[idx]); isDyn_final.append(isDyn2[idx]); isHqe_final.append(isHqe2[idx]); tts_final.append(tts2[idx])
            amp_final.append(amp2[idx]); hv_final.append(hv2[idx]); pvsv_final.append(pvsv2[idx]); svsn_final.append(svsn2[idx])
            riset_final.append(riset2[idx]); fallt_final.append(fallt2[idx]); fwhm_final.append(fwhm2[idx]); tts_ss_final.append(tts_ss_mcp[idx])
            num += 1
        elif i in hqmcp_id:
            idx = hqmcp_id.index(i)
            idx_final.append(hqmcp_id[idx])
            sn_final.append(sn3[idx]); gain_final.append(gain3[idx]); rsl_final.append(rsl3[idx]); dcr_final.append(dcr3[idx])
            pde_final.append(pde3[idx]); isDyn_final.append(isDyn3[idx]); isHqe_final.append(isHqe3[idx]); tts_final.append(tts3[idx])
            amp_final.append(amp3[idx]); hv_final.append(hv3[idx]); pvsv_final.append(pvsv3[idx]); svsn_final.append(svsn3[idx])
            riset_final.append(riset3[idx]); fallt_final.append(fallt3[idx]); fwhm_final.append(fwhm3[idx]); tts_ss_final.append(tts_ss_mcp[idx])
            num += 1
        else:
            print("==> No such PMT Index !")
            continue


    datadict = {"index": idx_final, "SN": sn_final, "isDyn": isDyn_final, "isHqe": isHqe_final, 
                "gain": gain_final, "resolution": rsl_final, "dcr": dcr_final, "pde": pde_final,
                "tts": tts_final, "amplitude": amp_final, "HV": hv_final, "PvsV": pvsv_final,
                "SvsN": svsn_final, "risetime": riset_final, "falltime": fallt_final, "fwhm": fwhm_final, 
                "tts_ss": tts_ss_final, "timeoffset": toff_final, "prePulseProb": ppp_final, 
                "afterPulseProb": app_final, "pmtx": pmtx_final, "pmty": pmty_final, "pmtz":pmtz_final}

    #datadict = {"index": idx_final, "SN": sn_final, "isDyn": isDyn_final, "isHqe": isHqe_final}
    datadf = pd.DataFrame(datadict)
    datadf.to_csv("PmtData_copy.csv")

    """

    
    read_pmtPos()

    """

