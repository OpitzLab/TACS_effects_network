# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:19:48 2023

@author: zhaoz
"""
def NetworkWeight(nettype):
    if nettype == 1: #alpha oscillation (10Hz)
        NetConnWeight = {
            "E2E_AMPA": 0.0012,
            "E2E_NMDA": 0.00003,
            "E2I_AMPA": 0.0004,
            "E2I_NMDA": 0.00003,
            "I2E_GABA": 0.0015,
            "I2I_GABA": 0.0004,
            "PoissonE": 0.00115,
            "PoissonI": 0.00006 
            }
    elif nettype == 2: #low beta oscillation (~15 Hz)
        NetConnWeight = {
            "E2E_AMPA": 0.0011,
            "E2E_NMDA": 0.00005,
            "E2I_AMPA": 0.0004,
            "E2I_NMDA": 0.00003,
            "I2E_GABA": 0.0015,
            "I2I_GABA": 0.0004,
            "PoissonE": 0.0010,
            "PoissonI": 0.00006    
            }
    elif nettype == 0: #no clear oscillation 
        NetConnWeight = {
            "E2E_AMPA": 0.0011,
            "E2E_NMDA": 0.00005,
            "E2I_AMPA": 0.0005,
            "E2I_NMDA": 0.00003,
            "I2E_GABA": 0.0015,
            "I2I_GABA": 0.0004,
            "PoissonE": 0.0012,
            "PoissonI": 0.00007 
            }
    return NetConnWeight