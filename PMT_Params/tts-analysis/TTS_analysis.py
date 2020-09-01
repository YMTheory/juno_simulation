#!/usr/bin/env python
# coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
from ROOT import TFile, TH1D, TH2D, TProfile, TLegend, TCanvas, TGraph, TGraphErrors
import uproot as up

import os
import sys

def getFiles(dir, suffix):
    """ collect all TTS scanning txt files """
    res = []
    for root, directory, files in os.walk(dir):
        for filename in files:
            name, suf = os.path.splitext(filename)
            if suf == suffix:
                res.append(os.path.join(root, filename))
    return res



def read_match(filename):
    """ read SN and scan_id corresponding and with PASS condition """
    data = pd.read_csv(filename, sep=" ")
    return np.array(data['Scan_ID'].loc[data["Indication_in_SS"]=="Pass"]), np.array(data['SN'].loc[data["Indication_in_SS"]=="Pass"])


def single_tube(argv):
    allFiles = getFiles("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/", ".txt")

    # read matches
    nnvt_id, nnvt_sn = read_match("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/instruction/valid_NNVT.csv")
    ham_id, ham_sn = read_match("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/instruction/valid_HAM.csv")


    graph_dyn = TGraphErrors()
    graph_mcp = TGraphErrors()

    nnvt_num, ham_num = 0, 0

    """ loop all files and generate pdf output """
    for idx, afile in enumerate(allFiles[0:100]):
        print("Processing scan file No {:d} -> {} " .format(idx, afile))
        #data1 = pd.read_csv(afile, header=[0], nrows=24, sep="\s+")  ## TTS
        data1 = pd.read_csv(afile, header=[25], nrows=24, sep='\s+') ## TT

        graph_dyn.Set(0)
        graph_mcp.Set(0)

        """ seperate MCP and Dynode PMTs """
        #scanid = int(re.sub("\D", "", afile))
        pattern = re.compile(r'(?<=scan)\d+')
        scanid = int(pattern.findall(afile)[0])

        if scanid in nnvt_id:
            print("It's a MCP PMT !!!")
            nnvt_num += 1
            
            """ average phi direction """
            #for itheta in range(1, 8):
            #    mean = np.array(data1.iloc[::, itheta]).mean()
            #    sigma = np.std( np.array(data1.iloc[::, itheta]) )
            #    graph_mcp.SetPoint(itheta-1, itheta, mean)
            #    graph_mcp.SetPointError(itheta-1, 0, sigma)
            
            """average theta direction"""
            for iphi in range(24):
                mean = np.array(data1.iloc[iphi, 1:8]).mean()
                sigma = np.std(np.array(data1.iloc[iphi, 1:8]))
                graph_mcp.SetPoint(iphi, iphi*15, mean)
                graph_mcp.SetPointError(iphi, 0, sigma)
            #plt.plot(graph_mcp.GetX(), graph_mcp.GetY(), "o-", ms=3)

        elif scanid in ham_id:
            print("It's a Dynode PMT !!!")
            ham_num += 1
            """ average phi direction """
            #for itheta in range(1, 8):
            #    mean = np.array(data1.iloc[::, itheta]).mean()
            #    sigma = np.std( np.array(data1.iloc[::, itheta]) )
            #    print("It's a Dynode PMT !!!")
            #    graph_dyn.SetPoint(itheta-1, itheta, mean)
            #    graph_dyn.SetPointError(itheta-1, 0, sigma)
            """average theta direction"""
            for iphi in range(24):
                mean = np.array(data1.iloc[iphi, 1:8]).mean()
                sigma = np.std(np.array(data1.iloc[iphi, 1:8]))
                graph_dyn.SetPoint(iphi, iphi*15, mean)
                graph_dyn.SetPointError(iphi, 0, sigma)
            plt.plot(graph_dyn.GetX(), graph_dyn.GetY(), "o-", ms=3)

    plt.title("Dyn PMT average all theta")
    plt.xlabel("phi")
    plt.ylabel("hittime/ns")
    plt.savefig(argv[1])

    print(" Dynode PMT Number: %d" % ham_num)
    print(" MCP PMT Number: %d" % nnvt_num)






def write_root(argv):
    allFiles = getFiles("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/", ".txt")
    led_array = ["led%s" %idx for idx in range(1,8)]
    print( "Total scan TTS files: {:d}" .format(len(allFiles)))

    # read matches
    nnvt_id, nnvt_sn = read_match("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/instruction/valid_NNVT.csv")
    ham_id, ham_sn = read_match("/Users/yumiao/Documents/Works/Simulation/PMT_Paramters/202008Updates/TTS-caimei-data-update/instruction/valid_HAM.csv")
    

    nnvt_num, ham_num = 0, 0
    # fill graph
    #hist_dyn = [TH2D("dyn%d"%i, "", 24, -7.5, 360-7.5, 100, 0, 1000) for i in range(7) ]
    #hist_nnvt = [TH2D("nnvt%d"%i, "", 24, -7.5, 360-7.5, 100, 0, 1000)  for i in range(7) ]
    hist_dyn = [TH1D("dyn%d"%i, "", 200, 300, 500) for i in range(7) ]
    hist_nnvt = [TH1D("nnvt%d"%i, "", 200, 300, 500)  for i in range(7) ]
    hist_dyn_phi = [ TH1D("dynPhi%d"%i, "", 200, 300, 500) for i in range(24) ]
    hist_mcp_phi = [TH1D("nnvtPhi%d"%i, "", 200, 300, 500) for i in range(24) ]
    hist_mcp_all = TH1D("nnvtall","", 200, 0, 1000)
    hist_dyn_all = TH1D("dynall","", 200, 0, 1000)

    #prof_dyn = [TProfile("dyn_prof%d"%i, "", 24, -7.5, 352.5, 0, 30) for i in range(7) ]
    #prof_nnvt = [TProfile("nnvt_prof%d"%i, "", 24, -7.5, 352.5, 0, 30) for i in range(7) ]
    
    hist2d_dyn = TH2D("dyn2d", "", 10, 0, 10, 100, 0, 30) 
    hist2d_nnvt = TH2D("nnvt2d", "", 10, 0, 10, 100, 0, 30)  

    """ loop all files and generate pdf output """
    for idx, afile in enumerate(allFiles):
        print("Processing scan file No {:d} -> {} " .format(idx, afile))
        #data1 = pd.read_csv(afile, header=[0], nrows=24, sep="\s+")  ## TTS 
        data1 = pd.read_csv(afile, header=[25], nrows=24, sep='\s+') ## TT

        """ seperate MCP and Dynode PMTs """
        #scanid = int(re.sub("\D", "", afile))
        pattern = re.compile(r'(?<=scan)\d+')
        scanid = int(pattern.findall(afile)[0])
        if scanid in nnvt_id:
            print("It's a MCP PMT !!!")
            nnvt_num += 1
            for idd, angle in enumerate(np.array(data1["angle"])):
                for col in range(1, 8):
                    #if np.array(data1.iloc[idd, col]) < 200000:
                    hist_mcp_all.Fill(np.array(data1.iloc[idd, col]))
                    if np.array(data1.iloc[idd, col]) > 380 and np.array(data1.iloc[idd, col]) <460:
                        hist_nnvt[col-1].Fill(np.array(data1.iloc[idd, col]))
                        hist_mcp_phi[idd].Fill(np.array(data1.iloc[idd, col]))
                        hist2d_nnvt.Fill(col, np.array(data1.iloc[idd, col]))
                        #hist_nnvt[col-1].Fill(angle, np.array(data1.iloc[idd, col]) )
                    #prof_nnvt[col-1].Fill(angle, np.array(data1.iloc[idd, col]) )
        elif scanid in ham_id:
            ham_num += 1
            print("It's a Dynode PMT !!!")
            for idd, angle in enumerate(np.array(data1["angle"])):
                for col in range(1, 8):
                    hist_dyn_all.Fill(np.array(data1.iloc[idd, col]))
                    #if np.array(data1.iloc[idd, col]) < 200000:
                    if np.array(data1.iloc[idd, col]) > 360 and np.array(data1.iloc[idd, col]) < 440:  # cut for station2 hittime test ...
                        hist_dyn[col-1].Fill(np.array(data1.iloc[idd, col]))
                        hist2d_dyn.Fill(col, np.array(data1.iloc[idd, col]))
                        hist_dyn_phi[idd].Fill(np.array(data1.iloc[idd, col]))
                    #hist_dyn[col-1].Fill(angle, np.array(data1.iloc[idd, col]) )
                    #prof_dyn[col-1].Fill(angle, np.array(data1.iloc[idd, col]) )

    print(" Dynode PMT Number: %d" % ham_num)
    print(" MCP PMT Number: %d" % nnvt_num)
 

    out = TFile(argv[1], "recreate")
    hist2d_dyn.Write()
    hist2d_nnvt.Write()
    hist_mcp_all.Write()
    hist_dyn_all.Write()
    for i in range(7):
        hist_dyn[i].Write()
        hist_nnvt[i].Write()
        #prof_dyn[i].Write()
        #prof_nnvt[i].Write()
    for j in range(24):
        hist_dyn_phi[j].Write()
        hist_mcp_phi[j].Write()
    out.Close()





def re_analysis(argv):
    """ re-analysis root file output """
    infile = TFile(argv[1], "read")
    #g_dyn = TGraph(); g_dyn.SetName("dyn")
    #g_mcp = TGraph(); g_mcp.SetName("mcp")
    dyn_arr, mcp_arr, dyn_err, mcp_err  = [], [], [], []
    dyn_phi, mcp_phi, dyn_phi_err, mcp_phi_err = [], [], [], []
    for i in range(7):
        h_dyn = infile.Get("dyn%d"%i)
        h_mcp = infile.Get("nnvt%d"%i)
        mean_dyn = h_dyn.GetMean()
        mean_mcp = h_mcp.GetMean()
        dyn_arr.append(mean_dyn)
        mcp_arr.append(mean_mcp)
        dyn_err.append(h_dyn.GetStdDev())
        mcp_err.append(h_mcp.GetStdDev())
    for j in range(24):
        h_dyn = infile.Get("dynPhi%d"%j)
        h_mcp = infile.Get("nnvtPhi%d"%j)
        mean_dyn = h_dyn.GetMean()
        mean_mcp = h_mcp.GetMean()
        dyn_phi.append(mean_dyn)
        mcp_phi.append(mean_mcp)
        dyn_phi_err.append(h_dyn.GetStdDev())
        mcp_phi_err.append(h_mcp.GetStdDev())

    """
    plt.plot([i for i in range(1, 8)], np.array(dyn_arr)-dyn_arr[0], "o-", color="blue", label="dynode")
    plt.plot([i for i in range(1, 8)], np.array(mcp_arr)-dyn_arr[0], "s-", color="hotpink", label="mcp")
    #plt.errorbar([i for i in range(1, 8)], np.sqrt(np.array(dyn_arr)**2-2.10**2), yerr=dyn_err, fmt="o-", color="blue", label="dynode")
    #plt.errorbar([i for i in range(1, 8)], np.sqrt(np.array(mcp_arr)**2-2.10**2), yerr=mcp_err, fmt="s-", color="hotpink", label="mcp")
    plt.legend()
    plt.xlabel("led No.")
    plt.ylabel("relative hittime mean/ns")
    plt.savefig("hittime_theta.pdf")

    """
    plt.plot([i*15 for i in range(24)], np.sqrt(np.array(dyn_phi)**2)-dyn_phi[0], "o-", color="blue", label="dynode")
    plt.plot([i*15 for i in range(24)], np.sqrt(np.array(mcp_phi)**2)-dyn_phi[0], "s-", color="hotpink", label="mcp")
    plt.legend()
    plt.xlabel("phi")
    plt.ylabel("hittime mean/ns")
    plt.savefig("hittime_phi.pdf")


def interpolate(argv):
    import scipy.interpolate as spi

    infile = TFile(argv[1], "read")
    #g_dyn = TGraph(); g_dyn.SetName("dyn")
    #g_mcp = TGraph(); g_mcp.SetName("mcp")
    dyn_arr, mcp_arr, dyn_err, mcp_err  = [], [], [], []
    dyn_phi, mcp_phi, dyn_phi_err, mcp_phi_err = [], [], [], []
    for i in range(7):
        h_dyn = infile.Get("dyn%d"%i)
        h_mcp = infile.Get("nnvt%d"%i)
        mean_dyn = h_dyn.GetMean()
        mean_mcp = h_mcp.GetMean()
        dyn_arr.append(mean_dyn)
        mcp_arr.append(mean_mcp)
        dyn_err.append(h_dyn.GetStdDev())
        mcp_err.append(h_mcp.GetStdDev())
    mcp_arr.insert(0, mcp_arr[0])
    mcp_arr.append(mcp_arr[-1])
    
    oldx_mcp = [0, 14, 30, 42.5, 55, 67, 77.5, 85, 90]
    ipo3 = spi.splrep(oldx_mcp, mcp_arr, k=2)
    newx = [i for i in range(0, 90)]
    iy3 = spi.splev(newx, ipo3)
    """
    from ROOT import TF1, TGraphErrors
    gMCP_tts = TGraphErrors()
    index = 0
    for x, y in zip(oldx_mcp, mcp_arr):
        gMCP_tts.SetPoint(index, x, y); index+=1
    ff = TF1("f1", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 90)
    gMCP_tts.Fit(ff, "R")
    cc = TCanvas()
    gMCP_tts.Draw("APL")
    cc.SaveAs("Fit.pdf")
    """
    plt.plot(oldx_mcp, mcp_arr, "o-", ms=4, color="blue", label="data")
    #plt.plot(newx, iy3, "-", color="salmon", label="interpolation")
    plt.legend()
    plt.xlabel("theta")
    plt.ylabel("TTS")
    plt.show()


if __name__ == "__main__":
    """ Need 2 parameters: first for PMT type(0 for dyn, 1 for nnvt), second is the output file name """
    #write_pdf(sys.argv)
    #write_root(sys.argv)
    #re_analysis(sys.argv)
    #single_tube(sys.argv)
    interpolate(sys.argv)
    
        

