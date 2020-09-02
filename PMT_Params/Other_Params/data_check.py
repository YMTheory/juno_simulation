#!/usr/bin/env python
# -*- coding=utf8 -*-
"""
# Author: MiaoYu ---> miaoyu@ihep.ac.cn
# Created Time : Wed Sep  2 18:52:14 2020
# File Name: data_check.py
"""

""" check pmt parameters generation """

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#plt.style.use("ggplot")

if __name__ == "__main__":
    df = pd.read_csv("PmtData.csv")

    pmtid = df["index"].to_numpy()
    isdyn = df["isDyn"].to_numpy()
    ishqe = df["isHqe"].to_numpy()
    tts_ss = df["tts_ss"].to_numpy()
    tts_ss_dyn = df["tts_ss"][df["isDyn"]==1].to_numpy()
    tts_ss_hqmcp = df["tts_ss"][(df["isDyn"]==0) & (df["isHqe"]==1)]
    tts_ss_lqmcp = df["tts_ss"][(df["isDyn"]==0) & (df["isHqe"]==0)]
    #plt.hist(tts_ss, bins=100, range=(0,20), color="darkorange", alpha=0.6,  label="allpmt: %.2f"%tts_ss.mean())
    plt.hist(tts_ss_dyn, bins=100,  range=(0,20), color="tomato", label="dynode: %.2f"%tts_ss_dyn.mean())
    plt.hist(tts_ss_hqmcp, bins=100, range=(0,20), color="skyblue", alpha=0.6,  label="hqmcp: %.2f"%tts_ss_hqmcp.mean())
    plt.hist(tts_ss_lqmcp, bins=100, range=(0,20), color="orchid", alpha=0.6,  label="lqmcp: %.2f"%tts_ss_lqmcp.mean())
    plt.xlabel("tts_ss/ns")
    plt.legend()
    plt.savefig("0902tts_ss.pdf")
    #plt.show()
