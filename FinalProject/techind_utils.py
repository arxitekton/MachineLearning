import os
import sys
import pandas as pd
import numpy as np 
import math
import time


# Commodity Channel Index 
def CCI(data, ndays): 
    TP = data['10Y Bond']
    name = 'CCI_' + str(ndays)
    CCI = pd.Series((TP - TP.rolling(window=ndays,center=False).mean()) / (0.015 * TP.rolling(window=ndays,center=False).std()), name = name)
    data = data.join(CCI)

    return data
 
    
# Moving average (MA)

# Simple Moving Average 
def SMA(data, ndays): 
    name = 'SMA_' + str(ndays) 
    SMA = pd.Series(data['10Y Bond'].rolling(window=ndays,center=False).mean(), name = name) 
    data = data.join(SMA) 
    
    return data


# Exponentially-weighted Moving Average 
def EWMA(data, ndays): 
    name = 'EWMA_' + str(ndays) 
    EMA = pd.Series(data['10Y Bond'].ewm(min_periods=(ndays - 1),span=ndays,ignore_na=False,adjust=True).mean(), name = name)    
    data = data.join(EMA)
    
    return data


# Rate of Change (ROC)
def ROC(data, ndays):
    name = 'ROC_' + str(ndays) 
    N = data['10Y Bond'].diff(ndays)
    D = data['10Y Bond'].shift(ndays)
    ROC = pd.Series(N/D, name=name)
    data = data.join(ROC)
    
    return data


# Compute the Bollinger Bands 
def BBANDS(data, ndays):
    uname = 'UpperBB_' + str(ndays)
    lname = 'LowerBB_' + str(ndays)
    MA = data['10Y Bond'].rolling(window=ndays).mean()
    SD = data['10Y Bond'].rolling(window=ndays).std()
    data[uname] = MA + (2 * SD) 
    data[lname] = MA - (2 * SD)
    
    return data