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

    toff = df["timeOffset"].to_numpy()
    toff_dyn = df["timeOffset"][df["MCP_Hama"]==1].to_numpy()
    toff_hqmcp = df["timeOffset"][(df["MCP_Hama"]==0) & (df["HiQE_MCP"]==1)]
    toff_lqmcp = df["timeOffset"][(df["MCP_Hama"]==0) & (df["HiQE_MCP"]==0)]
    #plt.hist(toff, bins=100, range=(0,100), color="darkorange", alpha=0.6,  label="allpmt: %.2f"%toff.mean())
    plt.hist(toff_dyn, bins=100,  range=(0,100), color="tomato", alpha=0.6, label="dynode: %.2f"%toff_dyn.mean())
    plt.hist(toff_hqmcp, bins=100, range=(0,100), color="skyblue", alpha=0.6,  label="hqmcp: %.2f"%toff_hqmcp.mean())
    plt.hist(toff_lqmcp, bins=100, range=(0,100), color="orchid", alpha=0.6,  label="lqmcp: %.2f"%toff_lqmcp.mean())
    plt.xlabel("toff")
    plt.legend()
    plt.savefig("0904toff.pdf")
    #plt.show()
