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
    df = pd.read_csv("PmtData_new.csv")

    dcr = df["DCR"].to_numpy()
    dcr_dyn = df["DCR"][df["MCP_Hama"]==1].to_numpy()
    dcr_hqmcp = df["DCR"][(df["MCP_Hama"]==0) & (df["HiQE_MCP"]==1)]
    dcr_lqmcp = df["DCR"][(df["MCP_Hama"]==0) & (df["HiQE_MCP"]==0)]
    plt.hist(dcr, bins=100, range=(0,100), color="darkorange", alpha=0.6,  label="allpmt: %.2f"%dcr.mean())
    plt.hist(dcr_dyn, bins=100,  range=(0,100), color="tomato", alpha=0.6, label="dynode: %.2f"%dcr_dyn.mean())
    plt.hist(dcr_hqmcp, bins=100, range=(0,100), color="skyblue", alpha=0.6,  label="hqmcp: %.2f"%dcr_hqmcp.mean())
    plt.hist(dcr_lqmcp, bins=100, range=(0,100), color="orchid", alpha=0.6,  label="lqmcp: %.2f"%dcr_lqmcp.mean())
    plt.xlabel("dcr")
    plt.legend()
    #plt.savefig("0904dcr.pdf")
    plt.show()
