# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 11:53:06 2020

@author: LilyHeAsamiko
"""
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as sci
import researchpy as rp
import pandas as pd
import statsmodels.formula.api as sfa
import statsmodels.api as sa
from sklearn.linear_model import LinearRegression
import quandl
import random
import scipy
import mpl_toolkits
from matplotlib import cm

def E(C_in, C_out, T, z):
#    R = 8.3144626 #J/K/mol
    R = 1.989 #Cal/K/mol
#    F = 96485.3329*10e-9 #C/mo
    F = 23.061 #kcal/V*gram
    E = R*T*np.log(C_out/C_in)/(F*10**3)
    return E

def I_CaL(V, GCaL):
    d = 1/(1+np.exp(-(V+7.24)/4.23))
    taud = 0.6+1/(np.exp(-0.05*(V+6)+np.exp(-0.09*(V+14))))
    finf = 1/(1+np.exp((V+22.88)/3.696))
    tauf = 7 + 1/(0.0045*np.exp(-(V+35.19)/10)+0.0045*np.exp(V+35.19/10))
    taus = 1000+1/(0.000035*np.exp(-(V+20.19)/4))+0.000035*np.exp(V+20.19/4)
    Aff = 0.6
    Afs = 1-Aff
    ff = 1
    fs = 1
    fCaf = 1
    fCas = 1
    f = Aff*ff+Afs*fs
    fCainf = finf
    taufCaf = 7 + 1/(0.004*np.exp(-(V-11.19)/7)+0.04*np.exp(V-11.19/7))
    taufCas = 100+1/(0.00012*np.exp(-(V+20)/3))+0.00012*np.exp(V+20/7)
    AfCaf = 0.3+0.6/(1+np.exp((V-10)/10))
    AfCas = 1-AfCaf
    fCa = AfCaf*fCaf+AfCas*fCas
    jCa = fCa
    tau = 75
    fCaMKinf = f
    fCaMK = fCaMKinf 

    taufCaMkf = 2.5*tauf
    AfCaMkf = Aff
    AfCaMks = 1-Aff
    fCaMKs = fs
    fCaMKf = ff
    fCaMk = AfCaMkf*fCaMKf+AfCaMks*fCaMKs
    fCaCaMkinf = f
    
    taufCaCaMkf = 2.5*taufCaMkf
    AfCaCaMkf = AfCaMkf
    AfCaCaMks = AfCaMks
    fCas = 1
    fCaCaMks = fCas
    fCaCaMk = AfCaCaMkf*fCaCaMkinf+AfCaCaMks*fCaCaMks
    Kmn = 0.002
    K2n = 1000
    K_2n = jCa
    cCasl = 1e-4
    alpha = 1/(K2n/K_2n+(1+np.exp(Kmn/cCasl))**4)
    gammaCai = 1
    gammaCao = 0.341
    zCa = 2
    VF = 86.55*23.0661
    R = 1.989 #Cal/K/mol
#    F = 96485.3329*10e-9 #C/mo
#    F = 23.061 #kcal/V*gram
    T = 300.15
    zCa = 2
    zK = 1
    zNa = 1
    cCass =  1e-4
    cCao = 1.8 #mM
#    GCaL = 7.7677e-5
    cNass = 8.23
    cNao = 140 #mM                                                                                                                                                                                                                                                                                                                                         bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
    cKss = 144
    cKo = 5.4 #mM
    PCa = 0.0001 #cm/s
    phiCa = zCa**2*VF**2/R/T*(gammaCai*cCass*np.exp(zCa*VF/R/T)-gammaCao*cCao)/(np.exp(zCa*VF/R/T)-1)
    ICaL_ = GCaL*phiCa
    PCaNa = 0.00125*PCa
    gammaNai = 0.75
    gammaNao = 0.75
    phiCaNa = zNa**2*VF**2/R/T*(gammaNai*cNass*np.exp(zNa*VF/R/T)-gammaNao*cNao)/(np.exp(zNa*VF/R/T)-1)
    ICaNa_ = PCaNa*phiCaNa
    PCaK = 3.574*10**(-4)*PCa
    gammaKi = 0.75
    gammaKo = 0.75
   # zK = 1
    phiCaK = zK**2*VF**2/R/T*(gammaKi*cKss*np.exp(zK*VF/R/T)-gammaKo*cKo)/(np.exp(zK*VF/R/T)-1)
    ICaK_ = PCaK*phiCaK
    PCaCaMK = 1.1*GCaL
    ICaLCaMK_ = PCaCaMK*phiCa
    PCaNaCaMK = 0.00125*PCaCaMK
    ICaNaCaMK_ = PCaNaCaMK*phiCaK   
    PCaKCaMK = 3.574*10**(-4)*PCaCaMK
    ICaKCaMK_ = PCaKCaMK*phiCaK
    KmCaMK = 0.15
    CaMK =  0.0124065 #sustained CaMKII is activated by CaMtrapping
    phiICaLCaMK = 1/(1+KmCaMK/CaMK)
    n = 0
    fCaMK = 1
    fCaCaMK = 1
    ICaL = ICaL_*d*(1-phiICaLCaMK)*(f*(1-n)+fCa*n*jCa)+ICaLCaMK_*d*phiICaLCaMK*(fCaMK*(1-n)*+fCaCaMK*n*jCa)
    ICaNa = ICaNa_*d*(1-phiICaLCaMK)*(f*(1-n)+fCa*n*jCa)+ICaNaCaMK_*d*phiICaLCaMK*(fCaMK*(1-n)*+fCaCaMK*n*jCa)
    ICaK = ICaK_*d*(1-phiICaLCaMK)*(f*(1-n)+fCa*n*jCa)+ICaKCaMK_*d*phiICaLCaMK*(fCaMK*(1-n)*+fCaCaMK*n*jCa)
    return [ICaL, ICaNa, ICaK]  

def I_Na(V, h_f, h_s,  G_Na, E_Na):
    minf = 1/(1+np.exp(-(V+48.6264)/9.871))
    taum = 1/(6.755*np.exp((V+11.64)/34.77)+8.552*np.exp(-(V+77.42)/5.955))
    hinf = 1/(1+np.exp((V+78.5)/6.22))
    tauh_f = 1/(3.686*10e-6*np.exp(-(V+3.8875)/7.8579)+16*np.exp((V-0.4963)/9.1843))
    tauh_s = 1/(0.009764*np.exp(-(V+17.95)/28.05)+0.3343*np.exp((V+5.730)/56.66))    
    Ah_f = 0.99
    Ah_s = 1 - Ah_f
    hCaMK_f = Ah_f*h_f
    hCaMK_s =  Ah_s*h_s
    jinf = hinf
    tauj = 4.8590 + 1/(0.8628*np.exp(-(V+166.7)/7.6)+1.1096*np.exp(-(V+6.2719)/9.0358))
    hCaMK_inf = 1/(1 + np.exp(-(V+84.7)/6.22))
    tauh_CaMK_s = 3*tauh_s
    Ah_CaMK_f = Ah_f
    Ah_CaMK_s = 1 - Ah_CaMK_f 
    h_CaMK_inf = h_f
    h_CaMK = Ah_CaMK_f*hCaMK_f + Ah_CaMK_s*hCaMK_s
    j_CaMK = jinf
    tauj_CaMK = 1.46*tauj
    Km_CaMK = 0.15
    phiINa_CaMK = 1/(1+Km_CaMK/(1+Km_CaMK))
    m = minf
    INa_f = G_Na*(V-E_Na)*m**3*((1-phiINa_CaMK)*h_CaMK*j_CaMK + phiINa_CaMK*h_CaMK*j_CaMK)
    return INa_f

def I_to(V, a, i_s, i_f, g, gto, EK):#transient Outward Current
    ainf = 1/np.exp((20-V)/13)
    taua = 1.0515/(1/(1.2*(1+np.exp(-(V-18.41)/29.38)))+3.5/(1+np.exp(-(V+100)/29.38)))
    i1_inf = 1/(1+np.exp((V+27)/13))
    i2_inf = i1_inf
    taui_S = 43 + 1/(0.001416*np.exp(-(V+96.52)/59.05)+1.7*10e-8*np.exp((V+11.41)/8.079))
    taui_F = 6.16 + 1/(0.39*np.exp(-(V+100)/100)+0.08*np.exp((V-8)/8.59))
    I_to = gto*a*i_s*i_f*g*(V-EK)
    return I_to

def I_sus(V, gsus, EK):
    a_sus = 1/(1+np.exp(-(V-12)/16))
    I_sus = gsus*a_sus*(V-EK)
    return I_sus

def I_Kr(V, CKpout, xr_s, xr_f, gkr, EK):
    xrinf = 1/(1+np.exp(-(V+8.337)/6.789))
    tauxr_f = 12.98 + 1/(0.3652*np.exp((V-14.06)/3.869)+4.123*10e-5*np.exp(-(V-30.18)/20.38))
    tauxr_s = 1.865 + 1/(0.06629*np.exp((V-19.7)/7.355)+1.128*10e-5*np.exp(-(V-12.54)/25.94))
    Axr_f = 1/(1+np.exp((V+54.81)/38.21))
    Axr_s = 1 - Axr_f 
    x_r = Axr_f*xr_f + Axr_s*xr_s 
    R_Kr = 1/(1+np.exp((V+55)/24))/(1+np.exp((V-10)/9.6))
    I_Kr = gkr*(CKpout/5.4)**0.5*x_r*R_Kr*(V-EK) 
    return I_Kr

def I_Ks(V, gks, Ca2p_sl, EK):
    xs1_inf = 1/(1+np.exp(-(V+11.6)/8.932))
    taux_s1 = 817.3 + 1/(0.001292*np.exp(-(V+210)/230)+2.326*10e-4*np.exp((V+48.28)/17.8))
    xs2_inf = xs1_inf   
    taux_s2 = 1/(0.01*np.exp((V-50)/20)+0.0193*np.exp(-(V+66.54)/31))
    I_Ks = gks*(1+0.6/(1+3.8*10e-5/Ca2p_sl)**1.4)*xs1_inf*xs2_inf*(V-EK) 
    return I_Ks

def I_f(V, gf_K, gf_Na, EK, ENa):
    y_inf = 1/(1+np.exp((V+87)/9.5))
    tauy = 2000/(np.exp(-(V+132)/10)+np.exp((V+57)/60))
    if_K = gf_K*y_inf*(V-EK)   
    if_Na = gf_Na*y_inf*(V-ENa) 
    I_f = if_K+if_Na
    return I_f

def I_K1(V, gk1, CKpout, Ek):
    xk1_inf = 1/(1+np.exp(-(V+2.5538*CKpout+144.59)/(1.5692*CKpout+3.8115)))
    taux_k1 = 122.2/(np.exp(-(V+127.2)/20.36)+np.exp((V+236.8)/69.33))
    R_K1 = 1/(1+np.exp((V+116-5.5*CKpout)/11))   
    taux_s2 = 1/(0.01*np.exp((V-50)/20)+0.0193*np.exp(-(V+66.54)/31))
    I_K1 = gk1*(CKpout/5.4)**0.5*xk1_inf*R_K1*(V-Ek) 
    return I_K1

def Calibration(V, fs):
    stair = np.repeat(0.1, 16760)
    stair[0:250] = np.repeat(-80, 250)
    stair[250:300] = np.linspace(-80, -120, 50)
    stair[300:700] = np.linspace(-120, -80, 400)
    stair[700:900] = np.repeat(-80, 200)
    stair[900:1900] = np.repeat(40, 1000)
    stair[1900:2400] = np.repeat(-120, 500)
    stair[2400:3400] = np.repeat(-80, 1000)
    
    stair[3400:3900] = np.repeat(-40, 500) 
    stair[3900:3920] = np.linspace(-40, -60, 20)
    stair[3920:4420] = np.repeat(-60, 500) 
    stair[4420:4440] = np.linspace(-60, -20, 20)
    stair[4440:4940] = np.repeat(-20, 500) 
    stair[4940:4960] = np.linspace(-20, -40, 20)
    stair[4960:5460] = np.repeat(-40, 500) 
    stair[5460:5480] = np.linspace(-40, 0, 20)
    stair[5480:5980] = np.repeat(0, 500) 
    stair[5980:6000] = np.linspace(0, -20, 20)
    stair[6000:6500] = np.repeat(-20, 500) 
    stair[6500:6520] = np.linspace(-20, 20, 20)
    stair[6520:7020] = np.repeat(20, 500) 
    stair[7020:7040] = np.linspace(20, 0, 20)
    stair[7040:7540] = np.repeat(0, 500) 
    stair[7540:7560] = np.linspace(0, 40, 20)
    stair[7560:8060] = np.repeat(40, 500) 
    stair[8060:8080] = np.linspace(40, 20, 20)
    stair[8080:8580] = np.repeat(20, 500) 
    stair[8580:14760] = list(stair[2400: 8580])[::-1]
    
    stair[14760:15260] = np.repeat(-80, 500)
    stair[15260:15760] = np.linspace(-80, 40, 500)
    stair[15760:15770] = np.linspace(40, -70, 10)
    stair[15770:15870] = np.linspace(-70, -110, 100)
    stair[15870:16260] = np.linspace(-110, -120, 390)
    stair[16260:16760] = np.repeat(-80, 500) 
    xnew = np.linspace(0,len(stair), round(len(stair)*fs/1000))
    xnew = np.linspace(0,len(stair), len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(stair)-1, len(stair)), stair)
#    temp = f(xnew)
    calib = f(xnew)

    return calib

def Val1(V):
    V1 = np.ones((1700, 7))+0.1
    t = 0
    while t < 7:
        V1[0:100, t] = np.repeat(-80, 100)
        V1[100: 1100, t] = np.repeat(-50+t*15, 1000)    
        V1[1100: 1600, t] = np.repeat(-40, 500)
        V1[1600: 1700, t] = np.repeat(-80, 100)
        t += 1 
    xnew = np.linspace(0,len(V1), len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(V1)-1, len(V1)), V1)
    V1 = f(xnew) 
    return V1

def Val2(V):
    V2 = np.ones((1200, 10))+0.1
    t = 0
    while t < 10:
        V2[0: 100, t] = np.repeat(-80, 100)
        V2[100: 600, t] = np.repeat(20, 500)        
        V2[600: 1100, t] = np.repeat(-140+t*20, 500)  
        V2[1100: 1200, t] = np.repeat(-80, 100)
        t += 1 
    xnew = np.linspace(0,len(V)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(V)-1, len(V)), V)
    V2 = f(xnew) 
    return V2

def Val3(V):
    V3 = np.repeat(0.1, 1350)
    V3[0: 100] = np.repeat(-80, 100)
    V3[100: 150] = np.repeat(-40, 50)        
    V3[150: 650] = np.repeat(20, 500)    
    V3[650: 1150] = np.repeat(-40, 500)
    V3[1150: 1350] = np.repeat(-80, 200)
    xnew = np.linspace(0,len(V)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(V), len(V)), V)
    V3 = f(xnew) 
    return V3

def Val4(V):
    DAD = np.repeat(0.1, 749)
    DAD[0: 50] = np.repeat(-80, 50)
    DAD[50: 53] = np.repeat(34, 3)        
    DAD[53: 61] = np.linspace(34, 30, 8)    
    DAD[61: 76] = np.linspace(30, 26, 15)
    DAD[76: 219] = np.linspace(26, -5, 143)    
    DAD[219: 257] = np.linspace(-5, -21, 38)    
    DAD[257: 325] = np.linspace(-21, -70, 68)
    DAD[325: 327] = np.repeat(-20, 2)    
    DAD[327: 347] = np.linspace(-20, -30, 20)    
    DAD[347: 357] = np.linspace(-30, -40, 10)
    DAD[357: 372] = np.linspace(-40, -65, 15)
    DAD[372: 384] = np.linspace(-65, -80, 12)
    DAD[384: 749] = np.repeat(-80, 365)
    xnew = np.linspace(0,len(DAD)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(DAD), len(DAD)), DAD)
    DAD_con = f(xnew) 
    return DAD_con

def Val5(V):
    EAD = np.repeat(0.1, 720)
    EAD[0: 50] = np.repeat(-80, 50)
    EAD[50: 53] = np.repeat(40, 3) 
    EAD[53: 56] = np.repeat(20, 3)       
    EAD[56: 76] = np.linspace(20, 30, 20)    
    EAD[76: 86] = np.repeat(30, 10)    
    EAD[86: 254] = np.linspace(30, -10, 168)    
    EAD[254: 305] = np.linspace(-10, -15.5, 51)
    EAD[305: 366] = np.linspace(-15.5, -20, 61)
    EAD[366: 426] = np.repeat(-20, 60)    
    EAD[426: 436] = np.linspace(-20, -10, 10)
    EAD[436: 446] = np.repeat(-10, 10)    
    EAD[446: 496] = np.linspace(-10, -20, 50)
    EAD[496: 516] = np.linspace(-20, -30, 20)
    EAD[516: 556] = np.linspace(-30, -75, 40)
    EAD[556: 606] = np.linspace(-75, -80, 50)
    EAD[606: 720] = np.repeat(-80, 114)
    xnew = np.linspace(0,len(EAD)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(EAD), len(EAD)), EAD)
    EAD_con = f(xnew) 
    return EAD_con

def Val6(V):
    AP = np.repeat(0.1, 395)
    AP[0: 50] = np.repeat(-80, 50)
    AP[50: 53] = np.repeat(34, 3) 
    AP[53: 61] = np.linspace(34, 30, 8)       
    AP[61: 76] = np.linspace(30, 26, 15)    
    AP[76: 260] = np.linspace(26, -8, 184)    
    AP[260: 299] = np.linspace(-8, -21, 39)    
    AP[299: 370] = np.linspace(-21, -68, 71)
    AP[370: 395] = np.linspace(-68, -80, 25)
    xnew = np.linspace(0,len(AP)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(AP), len(AP)), AP)
    AP_con = f(xnew) 
    return AP_con

def Val6(V):
    AP = np.repeat(0.1, 395)
    AP[0: 50] = np.repeat(-80, 50)
    AP[50: 53] = np.repeat(34, 3) 
    AP[53: 61] = np.linspace(34, 30, 8)       
    AP[61: 76] = np.linspace(30, 26, 15)    
    AP[76: 260] = np.linspace(26, -8, 184)    
    AP[260: 299] = np.linspace(-8, -21, 39)    
    AP[299: 370] = np.linspace(-21, -68, 71)
    AP[370: 395] = np.linspace(-68, -80, 25)
    xnew = np.linspace(0,len(AP)-1, len(V))
    f = scipy.interpolate.interp1d(np.linspace(0,len(AP), len(AP)), AP)
    AP_con = f(xnew) 
    return AP_con

def Step_ramp(V): #repeated 5s intervals
    SR = np.repeat(0.1, 10000)
    SR[0: 1000] = np.repeat(-80, 1000) # depolarize
    SR[1000: 3000] = np.repeat(20, 2000) # activate
    SR[3000: 5000] = np.repeat(-50, 2000) # repolarize
    SR[5000:10000] = np.repeat(-80, 5000) # deactivate
    return SR

def Step_pulse(V): #repeated 10s 
    SP = np.repeat(0.1, 5000)
    SP[0: 5000] = np.repeat(0.1, 5000)    
    SP[0: 1000] = np.repeat(-80, 1000)    
    SP[1000: 2000] = np.repeat(20, 1000)    
    SP[2000: 5000] = np.repeat(-80, 3000)
    return SP

def Extracellular_Unipolar_Potential(cpin, cpout, Vm, d, r):
    phi_e = r**2*cpin/4*cpout*sum(np.dot(-np.gradient(Vm),np.gradient(1/d)))
    return phi_e

def FV(D, Vm_x, X, I, dx):
    dVdt = D*(Vm_x[0:len(Vm_x)-2]-2*Vm_x[0:len(Vm_x)-1]+Vm_x[0:len(Vm_x)])/dx**2-I
    return dVdt

def Vpredict(Vt, Is, d, Ib, t):
# slice along n vector
    #   t = 10
#   Is = 1 muA
#    d = 0.5 
#    Ib = np.linspace(1.1, 2.1, 10)
#    Ib = range(1,3) transmembrane current
    x0 = 1.5*d
    sigma = 3
#    V = np.abs(V)
#    Vt = V[t]
    
    V_hat = []
    rhot_hat = []

    #1 lamdba = (rhom**b/(rhoi**b+rho0**b))**0.5 << d
    for b in range(len(Ib)):
        rhoi = abs(1/sigma/(0.75-0.85))#intracellular
        rho0 = abs(1/sigma/(1-0.75+0.85))#interstitial
#        rhot = (4*np.pi*d*V1[t]/Is/(1+rho0**Ib[b]/rhoi**Ib[b]*(2*np.exp(-d/Lambda)-np.exp(-2*d/Lambda))))**(1/Ib[b])
        rhot = (abs(4*np.pi*d*Vt/Is))**(1/Ib[b])
        Lambda = (rhot**Ib[b]/(abs(rhoi**Ib[b]+abs(rho0)**Ib[b])))**0.5
#        rho0_hat = (V1[t]*4*np.pi*d/Is/rhot**b*rhoi**b/(2*np.exp(-d/Lambda)-np.exp(-2*d/Lambda)))**1/b 
      #  if Labda >= d:
        V_hat.append(rhot**Ib[b]/4/np.pi/d*Is*(1+rho0**Ib[b]/rhoi**Ib[b]*(2*np.exp(-d/Lambda)-np.exp(-2*d/Lambda))))
        rhot_hat.append((4*np.pi*d*V_hat[-1]/Is)**(1/Ib[b]))
    #        rhoi_hat = 2*np.pi*d*V1[t]/Is
    #        rho0_hat = (V1[t]*4*np.pi*d/Is/rhot**b*rhoi**b/(2*np.exp(-d/Lambda)-np.exp(-2*d/Lambda)))**1/b
    #    else:
    #2 lambda >> d  , rho0 = 0
    #        rhoi = (4*np.pi*d*V1[t]/Is)**(1/b)        
    #        rhot = rhoi
    #        rho0 = 1/sigma/(1-0.75+0.85)
    #        Lambda = (rhot**b/(rhoi**b+rho0**b))**0.5    
    #        if Lambda > d:                    
        #        rhoi = 1/sigma/(0.75-0.85)
 
    #            V_hat = rhot_hat**b/4/np.pi/d*Is*(1+rho0_hat**b/rhoi**b*(2*np.exp(-d/Lambda)-np.exp(-2*d/Lambda)))
    return [V_hat, rhot_hat]
    
def ELBO_simplified(source0, dsource, V, Is, d, t,Ib,iters):
    #mesh 1.67*1.67*4, from source0 = (3.33,6.67,6.67, 0, 0, 0) to (60, 60, 60, ), dsource=(0.3, 1.9, 2.5, -3.7, 3.2, -34.7)
    #t = 10
    #z = 6.67
    #iters = 100
    V = np.array(V)
    source0 = [3.33,6.67,6.67, 0, 0, 0] 
    dsource = [0.3, 1.9, 2.5, -3.7, 3.2, -34.7]

    z = source0[2]
    x = source0[0]
    y = source0[1]
    alpha0 = source0[3]
    beta0 = source0[4] 
    theta = source0[5]
    dz = dsource[2]/180*np.pi
    dx = dsource[0]/180*np.pi
    dy = dsource[1]/180*np.pi
    dalpha = dsource[3]/100*alpha0
    dbeta = dsource[4]/100*alpha0
    dtheta = dsource[5]/100*alpha0
    V1_a = []
    P_a = []
    Q_a = []   
    temp= []
    while iters >0:
        temp = abs(V)
        p = []
        rho = []
        q = []
        for i in range(len(V)):
            #Add drift baseline b(t) = C*sum(a_k*cos(2*pi*df*t)+phi_k)
            temp[i] += random.gauss(0, 1)/len(temp)
        V1_a.append(temp)
        iters -= 1
        p = Vpredict(temp[t], Is, d, Ib, t)[0]
        rho = Vpredict(temp[t], Is, d, Ib, t)[1]
        q = temp[t]/np.cos(alpha0)*Is*np.array(rho) #Vi = Ji*cos(alpha0+dalpha)/Is/rho
        P_a.append(p)
        Q_a.append(q)
#    P_array = np.array(P_a)
#    Q_array = np.array(Q_a)    
#    P_param = np.mean(P_array[random.sample(range(100),20),:],0)  # q(theta|X), p(theta)  
    P_param =  np.mean(P_a,0)
    P_measure = P_a
    KLtemp = P_param*np.log(P_param-P_measure)
    KLtemp[np.isnan(KLtemp)]=0
    if np.size(np.shape(KLtemp))==2:
        #KL = sum(KLtemp, 1)
        #E = np.mean(P_measure,1)# E[P(X|theta)]
        Obj = np.mean(-KLtemp+P_measure,1)
        b = np.argmin(Obj)
        meanObj = []
        stdObj = []
        CV = []
        for i in range(1,len(Obj)):
            meanObj.append(np.mean(Obj[0:i]))
            stdObj.append(np.std(Obj[0:i]))
            CV.append(np.std(Obj[0:i])/np.mean(Obj[0:i]))
        plt.figure()
        plt.errorbar(np.linspace(0, len(Obj)-2, len(Obj)-1),meanObj, yerr = np.array(stdObj)/10)
        plt.title('timeseries_convergency: std')
        plt.figure()
        plt.errorbar(np.linspace(0, len(Obj)-2, len(Obj)-1),meanObj, yerr = np.array(CV)/10)
        plt.title('timeseries_convergency: CV')        
        Ib_opt = Ib[b]
        P = P_array[:, b]
        Q = Q_array[:, b]  
        plt.figure()
        plt.plot(P_a)
        plt.plot(Q_a)
        plt.title('reconstruction')
        plt.figure()
        plt.plot(abs(np.array(P_a)-np.array(Q_a)))
        plt.title('Residule')
        
        return [P, Q, Ib[b], Obj] 
    elif np.size(np.shape(KLtemp))==3:#iter, spatial, temporal 
        KL = np.mean(KLtemp, 0)
        E = np.mean(P_measure,0)# E[P(X|theta)]
        Obj = -KL+E
        b = np.argmin(np.sum(Obj,1))
        Ib_opt = Ib[b]
        P = P_array[:, b]
        Q = Q_array[:, b]  
        return  [P, Q, Ib, b, Obj[b]] 
    
    

def biomarkers(VECG):
    #Analysis: extract biomarkers: Tpeak, QT, rheart, Tamp, Tarea, TsymG, TsymA
    plt.figure()
    plt.plot(VECG)
    plt.title('One complete interval of ECG')
    tRpeak = np.argmax(VECG)
    Rpeak = np.max(VECG)
    tQon = np.argmin(VECG[0:tRpeak])
    tJ = np.argmin(VECG[tRpeak+10:])+tRpeak #83
    tTPeak = np.argmax(VECG[tRpeak+1:])+tQon #116
    TPeak = tTPeak - tJ #
    Tamp = max(VECG[tRpeak+1:])
    
    tToff = len(VECG)
    QT = tToff - tQon
    TsymG = abs(VECG[tTPeak+1]-[tJ])/abs(VECG[len(VECG)-1]-VECG[tTPeak+1])
    #TsymA = sum(VECG[range(tJ, tTPeak+1)])/sum(VECG[range(tTPeak+1,len(VECG))])
#    TsymG = (VECG[116]-VECG[114])/(VECG[117]-VECG[119])    
    TsymA = np.sum(VECG[tJ: tTPeak+1])/np.sum(VECG[tTPeak+1:])   
    return [tRpeak, Rpeak, tQon, tJ, tTPeak, TPeak, Tamp, tToff, QT, TsymA, TsymG]

def Validate_dynamicsVECG(VECG, VECG_new, samples):
    #tRpeak, Rpeak, tQon, tJ, tTPeak, TPeak, Tamp, tToff, QT, TsymA, TsymG
    origin = biomarkers(VECG)
    #tRpeak1, Rpeak1, tQon1, tJ1, tTPeak1, TPeak1, Tamp1, tToff1, QT1, TsymA1, TsymG1 = biomarkers(VECG_new)
    new = biomarkers(VECG_new)
    rel_dif=abs(np.array(new).ravel()-np.array(origin).ravel())/abs(np.array(origin).ravel())
#    rel_dif[np.isnan(rel_dif)] = 0 
    sum = 0
    for i in range(len(rel_dif)):
        if np.isnan(rel_dif[i]):
            rel_dif[i] = 0
        sum +=  np.mean(rel_dif[i])
    res = sum/np.shape(rel_dif)[0]
    if res <0.1:
        samples.append(VECG_new)
    return samples
            

def approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkr, V, EK, state):
    # one dimension V should be a verticle vector
    # 9 parameters units: s-1, V-1, s-1, V-1, ..., pS
    k1 = p1*np.exp(p2*V) #activation rate CI a*(1-r) I, C a*r O
    k2 = p3*np.exp(-p4*V) #deactivation rate I (1-r)*(1-a) CI, O r*(1-a) C
    k3 = p5*np.exp(p6*V) #inactivation rate O a*(1-r) I, C (1-r)*(1-a) CI
    k4 = p7*np.exp(-p8*V) #recovery rate I a*r O, CI r*(1-a) C
    ainf = k1/(k1 + k2)
    rinf = k4/(k3 + k4)
    taua = 1/(k1 + k2)
    taur = 1/(k3 + k4)
    da = ainf/taua*dt
    dr = rinf/taur*dt
#    da = ainf/taua - np.repeat(np.linspace(0, np.shape(V)[1]-1, np.shape(V)[1]),np.shape(V)[0]).reshape(np.shape(V))*(1/dt/taua +1/dt)    
#    dr = rinf/taur - np.repeat(np.linspace(0, np.shape(V)[1]-1, np.shape(V)[1]),np.shape(V)[0]).reshape(np.shape(V))*(1/dt/taur +1/dt)    
    if np.shape(V)[0] == 1:
#        print()
        exit('input one column of time series')
    elif np.size(V) == np.shape(V)[0]:  
        a = ainf -da*np.linspace(0, np.shape(V)[0]-1, np.shape(V)[0])/dt
        r = rinf -dr*np.linspace(0, np.shape(V)[0]-1, np.shape(V)[0])/dt
    else:      
        a = ainf - da*np.repeat(np.linspace(0, np.shape(V)[0]-1, np.shape(V)[0])/dt, np.shape(V)[1]).reshape(np.shape(V))
        r = rinf - dr*np.repeat(np.linspace(0, np.shape(V)[0]-1, np.shape(V)[0])/dt, np.shape(V)[1]).reshape(np.shape(V))

    if state == 1:
        Ikr = gkr*a*r*(V-EK)
    elif state == 2:
        Ikr = gkr*a*(1-r)*(V-EK) 
    elif state == 3:
        Ikr = gkr*(1-a)*(1-r)*(V-EK) 
    elif state == 4:
        Ikr = gkr*(1-a)*r*(V-EK)         
    return Ikr


def Calibrate_Validate(I0, calib, gkr0, tp, N, states, drifts, dt, p1, p2, p3, p4, p5, p6, p7, p8, EK):
    EXP = np.exp(I0)
    A = np.sort(np.linspace(min(I0), max(I0), N))
    B = np.sort(np.linspace(min(EXP), max(EXP), N)/max(EXP)/10)
    dVmax = 160
    p1 = 1/(450*dt)
#    p2 = 1/40/dVmax
    p2 = -1/40
    p3 = 1/(500*dt)
    p4 = -1/160
    p5 = 1/(500*dt)
    p6 = 1/120
    p7 = 1/(500*dt)
    p8 = -1/120
    IkrO = []
    GkrO = []
    ResiduleO = []
    Ikr = []
    Gkr = []
    Residule = []
    
    for i in range(len(tp)):
        V = calib[(tp[i]+drifts[i]):(tp[i]+drifts[i]+N)] 
        V1 = Val1(V)
        V2 = Val2(V)
        V3 = Val3(V)
        DAD = Val4(V)
        EAD = Val5(V)
        AP = Val6(V)
        state = states[i]
        I0temp = (I0[tp[i]:tp[i]+N]).ravel()
        gkr = gkr0
        Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkr, V, EK, state), dtype = float)
        gkr = I0temp/Itemp*gkr 
        Itempcorr = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkr, V, EK, state), dtype = float)

        residule = Itempcorr-I0temp
        gkrcor = I0temp/Itempcorr*gkr
        IkrO.append(Itempcorr)
        GkrO.append(gkrcor)
        ResiduleO.append(residule)
        
        Itempval = Itempcorr
        residuleval  = residule
        gkrval = gkrcor
        

        #Gkr.append(np.mean(gkrcor))
        #validation 1  
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, np.repeat(gkrcor[0:min(N, 1700)], np.shape(V1)[1]).reshape(min(N, 1700), np.shape(V1)[1]), -V1[0:min(N, 1700)], EK, state), dtype = float) #Val2(V[0:1700])
        else:  
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, np.repeat(gkrcor[0:min(N, 1700)], np.shape(V1)[1]).reshape(min(N, 1700), np.shape(V1)[1]), V1[0:min(N, 1700)], EK, state), dtype = float) #Val2(V[0:1700])
        gkrtemp = np.repeat(I0temp[0:min(N, np.shape(V1)[0])]/Itemp[:-1], np.shape(V1)[1]).reshape(min(N, np.shape(V1)[0]), np.shape(V1)[1])*np.repeat(gkr[0:min(N, np.shape(V1)[0])], np.shape(V1)[1]).reshape(min(N, np.shape(V1)[0]), np.shape(V1)[1])    
        
        if sum(abs(Itemp[:,-1]-I0temp[0:min(N, np.shape(V1)[0])]) <= residuleval[0:min(N, np.shape(V1)[0])])>0:
            idxtemp = abs(Itemp[:,-1]-I0temp[0:min(N, np.shape(V1)[0])]) <= residuleval[0:min(N, np.shape(V1)[0])]
            Itempval[idxtemp] = Itemp[idxtemp, -1]
            gkrval[idxtemp] = gkrtemp[idxtemp, -1] 
            residuleval = Itempval-I0temp
        
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation1')
    
            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')  
            plt.title(str(i)+' '+'th'+' '+'residule for validation1')
        #    if sum(residuleval <= residule)>=0:
        #        residule = residuleval
        #validation 2  
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, np.repeat(gkrcor[0:min(N, np.shape(V2)[0])], np.shape(V2)[1]).reshape(min(N, np.shape(V2)[0]), np.shape(V2)[1]), -V2[0:min(N, np.shape(V2)[0])], EK, state), dtype = float) #Val2(V[0:1200])
        else:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, np.repeat(gkrcor[0:min(N, np.shape(V2)[0])], np.shape(V2)[1]).reshape(min(N, np.shape(V2)[0]), np.shape(V2)[1]), V2[0:min(N, np.shape(V2)[0])], EK, state), dtype = float) #Val2(V[0:1200])
        gkrtemp = np.repeat(I0temp[0:min(N, np.shape(V2)[0])]/Itemp[:,-1], np.shape(V2)[1]).reshape(min(N, np.shape(V2)[0]), np.shape(V2)[1])*np.repeat(gkr[0:min(N, np.shape(V2)[0])], np.shape(V2)[1]).reshape(min(N, np.shape(V2)[0]), np.shape(V2)[1])    
        if sum(abs(Itemp[:,-1]-I0temp[0:min(N, np.shape(V2)[0])]) <= residuleval[0:min(N, np.shape(V2)[0])])>0:
            idxtemp = abs(Itemp[:,-1]-I0temp[0:min(N, np.shape(V2)[0])]) <= residuleval[0:min(N, np.shape(V2)[0])]
            for i in range(len(Itempval)):
                if i >= len(idxtemp):
                    break
                if idxtemp[i]== True:
                    Itempval[i,] = Itemp[i, -1]                
                    gkrval[i,] = gkrtemp[i, -1] 
            residuleval = Itempval-I0temp
            
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation2')

            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')
            plt.title(str(i)+' '+'th'+' '+'residule for validation2')
            
        #validation 3  
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(V3)[0])], -V3[0:min(N, np.shape(V3)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        else:    
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(V3)[0])], V3[0:min(N, np.shape(V3)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        gkrtemp = I0temp[0:min(N, np.shape(V3)[0])]/Itemp*gkr[0:min(N, np.shape(V3)[0])] 
        if sum(abs(Itemp-I0temp[0:min(N, np.shape(V3)[0])]) <= residuleval[0:min(N, np.shape(V3)[0])])>0:
            idxtemp = abs(Itemp-I0temp[0:min(N, np.shape(V3)[0])]) <= residuleval[0:min(N, np.shape(V3)[0])]
            for i in range(len(Itempval)):
                if i >= len(idxtemp):
                    break
                if idxtemp[i]== True:
                    Itempval[i] = Itemp[i]                
                    gkrval[i] = gkrtemp[i] 
            residuleval = Itempval-I0temp
            
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation3')

            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')
            plt.title(str(i)+' '+'th'+' '+'residule for validation 3')
            
        #validation 4  
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(DAD)[0])], -DAD[0:min(N, np.shape(DAD)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        else:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(DAD)[0])], DAD[0:min(N, np.shape(DAD)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        gkrtemp = I0temp[0:min(N, np.shape(DAD)[0])]/Itemp*gkr[0:min(N, np.shape(DAD)[0])] 
        if sum(abs(Itemp-I0temp[0:min(N, np.shape(DAD)[0])]) <= residuleval[0:min(N, np.shape(DAD)[0])])>0:
            idxtemp = abs(Itemp-I0temp[0:min(N, np.shape(DAD)[0])]) <= residuleval[0:min(N, np.shape(DAD)[0])]
            for i in range(len(Itempval)):
                if i >= len(idxtemp):
                    break
                if idxtemp[i]== True:
                    Itempval[i] = Itemp[i]                
                    gkrval[i] = gkrtemp[i] 
            residuleval = Itempval-I0temp
            
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation4')

            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')  
            plt.title(str(i)+' '+'th'+' '+'residule for validation4')
        
        #validation 5 
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(EAD)[0])], -EAD[0:min(N, np.shape(EAD)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        else:
            Itemp = abs(Itemp-I0temp[0:min(N, np.shape(EAD)[0])]) <= residuleval[0:min(N, np.shape(EAD)[0])] >0
            idxtemp = ap.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(EAD)[0])], EAD[0:min(N, np.shape(EAD)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        gkrtemp = I0temp[0:min(N, np.shape(EAD)[0])]/Itemp*gkr[0:min(N, np.shape(EAD)[0])] 
        if sum(abs(Ibs(Itemp-I0temp[0:min(N, np.shape(EAD)[0])]))) <= residuleval[0:min(N, np.shape(EAD)[0])]
            for i in range(len(Itempval)):
                if i >= len(idxtemp):
                    break
                if idxtemp[i]== True:
                    Itempval[i] = Itemp[i]                
                    gkrval[i] = gkrtemp[i] 
            residuleval = Itempval-I0temp
            
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation5')

            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')  
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation5')
            
        #validation 6  
        if min(I0temp) <0:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(AP)[0])], -AP[0:min(N, np.shape(AP)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        else:
            Itemp = np.array(approximate(dt, p1, p2, p3, p4, p5, p6, p7, p8, gkrcor[0:min(N, np.shape(AP)[0])], AP[0:min(N, np.shape(AP)[0])], EK, state), dtype = float) #Val2(V[0:1700])
        gkrtemp = I0temp[0:min(N, np.shape(AP)[0])]/Itemp*gkr[0:min(N, np.shape(AP)[0])] 
        if sum(abs(Itemp-I0temp[0:min(N, np.shape(AP)[0])]) <= residuleval[0:min(N, np.shape(AP)[0])])>0:
            idxtemp = abs(Itemp-I0temp[0:min(N, np.shape(AP)[0])]) <= residuleval[0:min(N, np.shape(AP)[0])]
            for i in range(len(Itempval)):
                if i >= len(idxtemp):
                    break
                if idxtemp[i]== True:
                    Itempval[i] = Itemp[i]                
                    gkrval[i] = gkrtemp[i] 
            residuleval = Itempval-I0temp
            
            plt.figure()
            plt.plot(Itempcorr,'r-')
            plt.plot(I0temp,'g--')
            plt.plot(Itempval, 'b*')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation6')

            plt.figure()
            plt.plot(residule,'r')
            plt.plot(residuleval, 'g')
            plt.title(str(i)+' '+'th'+' '+'I reconstruction for validation6')
        
        Ikr.append(Itempval)
        Gkr.append(gkrval)
        Residule.append(residuleval)
        #
    return [IkrO, GkrO, ResiduleO, Ikr, Gkr, Residule]        

def Sensitivity(B1, B2):#B1: initiate, #B2: scaled
    dBMax = np.max([np.max(B1, 1), np.max(B2, 1)], 0)
    dBMin = np.min([np.min(B1, 1), np.min(B2, 1)], 0)
    S = (B2 -B1)/np.max(abs(dBMax)-abs(dBMin))
    return S

def main():
    E0 = pd.read_csv('D:/ECG/herg25oc1-staircaseramp-A04.csv', header = [0]) 
    E1 = pd.read_csv('D:/ECG/herg25oc1-staircaseramp-A01-after.csv', header = [0])
    E0Arr = np.array(E0, dtype = float)
    plt.figure()
    plt.plot(E0Arr)
    E0f = np.fft.fft(E0Arr)
    E0f = np.mean(E0f)    
    E0f = np.fft.fft(E0Arr)-np.mean(E0f)
    plt.figure()    
    plt.plot(E0f)    
    E0_cor = np.fft.ifft(E0f[1:])
    plt.figure()
    plt.plot(E0_cor) 

    # AnotherTest withAp06: CdV/dt = -(ICaL+ItoX+ItoY+Ik1+IKs+Ikr+If+IK+INa+ICa+Iothers-Istim)
    # Data:
    # BayesInferenced(Default with Abstract/Normal Distribution)
    # Drug:
    # Period:
    # Duration:
    # QT interval(Toff -Qon): 
    # on basis of AP: Ic50 muM(raw: pIC50 M), Hill
    # Validation with DAD I, EAD II(500, 250)
    # Calibration with AP, or CaTr or AP and CaTr
    # source of cells in IonicChannals: Slow outward: Iks(Slow delayed rectifier potassium)+ ICaL()+ INa(Fast-Sodium) +  Ito(transient to outward), slow inward(Ik1, Inward Rectifier Potassium)+ fast inward: IKr(Rapid Delayed Rectifier Potassium)
    # Endocardial: 200 Ikr, Iks, If, ICaL
    # Mid-myocardial: 150 AP, ItoX
    # Epicardial: 125 Ik, INa, ICa, ItoY,Ik1, 
    # postprocess: AP rate, Peak, Upstroke phase#, ERD/LRD30, 60
    # Ikr HERG expressed (rapidly inactivating delayed rectifier)
    # AP BCL: 1000ms after 
    # Phase 0: resting
    # Phase 1: depolarisation
    # Phase 2: plateau
    # Phase 3: repolarisation
    # drug-induced prolong: APD increases the risk of early afterdepolarisation(EAD):secondary depolarisation
    # occurs at 2/3 (P28, 1-4), might critical: Tpeak-Toff
    # 1. initiate an ectopic beat 
    # 2. (Given right condition/substrate) trigger a wave of activation  
    # 3. degenerate into torsades 
    fs = 4680
    EAD = pd.read_csv('C:/Users/LilyHeAsamiko/VM Shared Files_bu/Chaste/heart/test/data/sample_APs/Ead.dat', header = [0]) 
    plt.figure()  
    plt.plot(EAD)
    plt.title('EAD_Data')   
    
    E0 = pd.read_csv('C:/Users/LilyHeAsamiko/VM Shared Files_bu/Chaste/heart/test/data/ionicmodels/FoxRegularStimValidData.csv', header = [0]) 
    AP = E0.iloc[:,1]
    Vm = AP
    plt.figure()  
    plt.plot(AP)
    plt.title('AP_Data')    
    Ik_dat = E0.iloc[:,2]
    plt.figure()  
    plt.plot(Ik_dat) 
    plt.title('AP_Data_m(K)')    
    INa_dat = E0.iloc[:,3]    
    plt.figure() 
    plt.plot(INa_dat)
    plt.title('AP_Data_h(Na)') 
    ICa_dat = E0.iloc[:,4]        
    plt.figure()  
    plt.plot(ICa_dat)
    plt.title('AP_Data_j(Ca)') 
    Ikr_dat = E0.iloc[:,5]      
    plt.figure()  
    plt.plot(Ikr_dat) 
    plt.title('AP_Data_kr')   
    Iks_dat = E0.iloc[:,6]     
    plt.figure()  
    plt.plot(Iks_dat)
    plt.title('AP_Data_ks') 
    ItoX_dat = E0.iloc[:,7]        
    plt.figure()  
    plt.plot(ItoX_dat)
    plt.title('AP_Data_toX')
    ItoY_dat = E0.iloc[:,8]         
    plt.figure()  
    plt.plot(ItoY_dat)
    plt.title('AP_Data_toY') 
    If_dat = E0.iloc[:,9]       
    plt.figure()  
    plt.plot(If_dat)   
    plt.title('AP_Data_f') 
    Ik1_dat = E0.iloc[:,10]     
    plt.figure()  
    plt.plot(Ik1_dat) 
    plt.title('AP_Data_k1')   
    ICaL_dat = E0.iloc[:,11]     
    plt.figure()  
    plt.plot(ICaL_dat)
    plt.title('AP_Data_CaL') 
#    Some of the Iothers:     
#    plt.figure()  
#    plt.plot(E0.iloc[:,12])
#    plt.title('AP_Data_CaI') 
#    plt.figure()  
#    plt.plot(E0.iloc[:,13])
#    plt.title('AP_Data_CaSr')    
       
    EE = pd.read_csv('C:/Users/LilyHeAsamiko/VM Shared Files_bu/Chaste/mesh/test/data/butterflyEle.csv', header = None,  engine = 'python') 
    E1 = pd.read_csv('C:/Users/LilyHeAsamiko/VM Shared Files_bu/Chaste/mesh/test/data/butterflyNode.csv', header = None)

#    IKr0_X = E1.iloc[min(np.array(EE.iloc[:,0], dtype = int),len(E1)), 0]*0.2
#    IKr0_Y = E1.iloc[min(np.array(EE.iloc[:,1], dtype = int),len(E2)), 1]*0.2    
    V0 = E1.iloc[:, 0]*0.2
    V1 = E1.iloc[:, 1]*0.2    
    
    plt.figure()  
    plt.plot(V0)
    plt.title('RawData_V0')    
    
    
    plt.figure()    
    plt.plot(V1)    
    plt.title('RawData_V1') 
    
    #preprocess on IKro_X
    E1Arr = np.array(V0[450:], dtype = float)
    plt.figure()
    plt.plot(E1Arr)
    
    E1Arr_avg = E1Arr
    for i in range(len(E1Arr)-20):
        E1Arr_avg[i] = np.mean(E1Arr[i:i+20])
    plt.plot(E1Arr_avg)
    plt.title('smoothed_V0')
                                                                                                  
    E1f = np.fft.fft(E1Arr_avg)
    E1f = np.mean(E1f)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    E1f = np.fft.fft(E1Arr_avg)-np.mean(E1f)

    plt.figure()    
    plt.plot(E1f)    
    E1_cor = np.fft.ifft(E1f[1:])
    plt.figure()
    plt.plot(E1_cor)  
    plt.title('Baseline Corrected preprocessed V0')    
    #endo: mid: epi : 5:3:2
    AP_endo = E1_cor[0:round(len(E1_cor)/10)*5]
    AP_mid = E1_cor[round(len(E1_cor)/10*5+1):round(len(E1_cor)/10*8)]
    AP_epi = E1_cor[round(len(E1_cor)/10*8+1):]

#    AP_endo_valid = []
#    AP_mid_valid = []
#    AP_epi_valid = []    
#    while i in range(len(AP_endo)):
#        temp = AP_endo[i+1:i+5000]
#        if abs(np.min(temp)-AP_endo[i]) <= 2/3*abs(np.max(temp)-AP_endo[i]):
#            print(i)
#            print(abs(np.min(temp)-AP_endo[i]))
#            print(abs(np.max(temp)-AP_endo[i]))
#            AP_endo_valid.append(AP_endo[i:np.argmin(temp)])
#            i = np.argmin(temp)+1
#        else:
#            i += 1
#    while i in range(len(AP_mid)):
#        temp = AP_mid[i+1:i+5000]
#        if abs(np.min(temp)-AP_mid[i]) <= abs(np.max(temp)-AP_mid[i])*2/3:
#            print(abs(np.min(temp)-AP_endo[i]))
#            print(abs(np.max(temp)-AP_endo[i]))
#            AP_mid_valid.append(AP_mid[i:np.argmin(temp)])            
#            i = np.argmin(temp)+1
#        else:
#            i += 1
#    while i in range(len(AP_epi)):
#        temp = AP_epi[i+1:i+5000]
#        if abs(np.min(temp)-AP_epi[i]) <= abs(np.max(temp)-AP_epi[i])*2/3:
#            print(abs(np.min(temp)-AP_endo[i]))
#            print(abs(np.max(temp)-AP_endo[i]))
#            AP_epi_valid.append(AP_epi[i:np.argmin(temp)])
#            i = np.argmin(temp)+1
#        else:
#            i += 1
    
       
    #check I
    a = 1
    i_f = 0.6
    i_s = 1-i_f
    xr_f = 0
    xr_s = 0.6
    g = 1
    gto = 0.192
    gkr = 0.0342
    gks = 0.0029
#    gsus = 0.0301
    gk1 = 0.04555
    gf_Na = 0.0116
    gf_K = 0.0232
    gCaL = 7.7677e-5
    fs = 4680 #Hz
    #Na block
    h_f = 0.8
    h_s = 0.8
    G_Na = 39.46
#    V = E0_cor[0:fs]
#    V = np.real(E1_cor[0:round(fs/2)])
#    V = np.real(AP_endo)
#    plt.plot(V)
#    plt.title('AP_endo_Data')    
    #IKr Block
    CKpin = 110 #mM
    CKpout = 5.4 #mM
    CCapin = 4.36e-5
    CCapout = 1.8    
    CNapin = 8.23
    CNapout= 140
    Ca2p_sl = 1e-4
    T = 300.15 #K
    

    EK = E(CKpin, CKpout, T, z=1)
    ECaL = E(CCapin, CCapout, T, z = 2) 
    ENa = E(CNapin, CNapout, T, z = 1)


    f = open("D:\PhD in Oxford,Ethz,KI,others\OxfordPhD\ECG\Butterfly_mesh_FHN\PsudoEcgGromElectrodeAt_-1_0_0.dat", "r")
    l = []
    for x in f:
        s = f.readline()
        s =s.split('\t')
        for string in s:
            l.append(string.split())
    f.close()
    
    lArr=np.array(l[:])
    n_ECG = int(len(lArr)/3)
    ECG = np.repeat(0.1, n_ECG)#filter_data
    ECG_r = np.repeat(0.1, n_ECG)#raw_data
    RR = 20/13 #s
    nRR = int(np.round(n_ECG/12-70))   
    #using zero baseline, filtered ECG
    for i in range(2, len(lArr), 3):
        ECG[int(i/3)] = float(lArr[i])  
        ECG_r[int(i/3)] = float(lArr[i-1]) 
    plt.figure()     
    plt.plot(ECG)  
    plt.plot(ECG_r)  
   
    
    Is = 1 #muA
    drho = 0.5 
    Ib = ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim
    #Ib = np.linspace(1.1, 2.1, 10)#transmembrane current

    t = 10
    VECG = ECG[188:188+nRR]
    ECG_Vt = VECG[t]
    
    Pred = Vpredict(ECG_Vt, Is, drho, Ib, t)
    VPred = Pred[0]
    RohPred = Pred[1]
    iters = 100
    
    source0 = [3.33, 6.67, 6.67, 0, 0, 0]
    dsource = [0.3, 1.9, 2.5, -3.7, 3.2, -34.7]
    #test on t = 10
    reconstructOut = ELBO_simplified(source0, dsource, V, Is, drho, t,Ib,iters)    
    P_rec = reconstructOut[0] 
    Q_rec = reconstructOut[1]
    Ib = reconstructOut[2] 
    b= reconstructOut[3] 
    KL = reconstructOut[4]
    res = []
    for i in range(1, iters):
        P_avg = np.mean(P_rec[0: i, ], 0)
        Q_avg = np.mean(Q_rec[0: i, ], 0)
        res.append(np.mean(abs(P_avg-Q_avg)))
    Iter_opt = np.argmin(res)
    P_rec = []
    plt.figure()
    plt.plot(res)
    plt.title('PQresidule_iters')
    plt.figure()
    plt.plot(P[Iter_opt,])
    plt.plot(Q[Iter_opt,])  
    plt.legend()    
    plt.title('Best reconstruction Q of P')   

'''
    #Analysis: extract biomarkers: Tpeak, QT, rheart, Tamp, Tarea, TsymG, TsymA
    plt.figure()
    plt.plot(VECG)
    plt.title('One complete interval of ECG')
    tRpeak = np.argmax(VECG)
    Rpeak = np.max(VECG)
    tQon = np.argmin(VECG[0:tRpeak])
    tJ = np.argmin(VECG[tRpeak+10:])+tRpeak #83
    tTPeak = np.argmax(VECG[tRpeak+1:])+tQon #116
    TPeak = tTPeak - tJ #
    Tamp = max(VECG[tRpeak+1:])
    
    tToff = len(VECG)
    QT = tToff - tQon
    #TsymA = sum(VECG[tJ:tTPeak+1])/sum(VECG[tTPeak+1:])
    #TsymA = sum(VECG[range(tJ, tTPeak+1)])/sum(VECG[range(tTPeak+1,len(VECG))])
    # TsymG = np.mean(VECG[84: 116]-VECG[83: 115])/np.mean(VECG[116:]-VECG[115:len(VECG)-1])
    TsymG = (VECG[116]-VECG[114])/(VECG[117]-VECG[119])    
    TsymA = np.sum(VECG[84: 116])/np.sum(VECG[116:])
'''    
    VECG_Samples = []
    iters = 100
    VECG_Samples.append(VECG)
    for i in range(iters):
        l = np.shape(VECG_Samples)[0]
        VECG_new = VECG + np.random.normal(0, 1, np.size(VECG))/100
        VECG_Samples = Validate_dynamicsVECG(VECG, VECG_new, samples)       
        if np.shape(VECG_Samples)[0] > l:
    random()

    print()
    d = np.sqrt(abs(IKr0_X-IKr0_Y)**2)
    r = 0.0175
    #phi_xy = Extracellular_Unipolar_Potential(CKpin, CKpout,  np.array([IKr0_X, IKr0_Y]), d, r)
#    Vm = np.array([IKr0_X, IKr0_Y])
    phi_e = r**2*CKpin/4*CKpout*sum(-np.dot(np.gradient(Vm),np.gradient(1/d)))
    plt.plot(phi_e)
    
    Ito = []
    ItoX_t = []
    ResItoX_t = []
    ItoY_t = []
    ResItoY_t = []
    I_Ks = []
    Iks_t = []
    ResIks_t = []
    If = []
    If_t = []
    ResIf_t = []
    Ik1 = []
    Ik1_t = []
    ResIk1_t = []
    ICaALL = []
    ICaALL_t = []
    ResICaALL_t = []
    INa_f = []
    INa_f_t = []
    ResINa_f_t = []
    Ikr = []
    Ikr_t = []
    ResIkr_t = []
    
    #VECG
    Iothers =  ICa_dat +E0.iloc[:,12]+E0.iloc[:,13]
    
    h = plt.hist(ItoX_dat, density = True, bins = 50)
    #AP = VECG
     AP = np.array(AP)
     EK = np.array(EK)
     ENa = np.array(ENa)
    # ECa = np.array(ECa)
    Ito = I_to(AP, a, i_s, i_f, g, gto, EK)
#    ItoX_t = I_to(AP, a, i_s, i_f, g, gto, EK)
    ItoX_t.append(I_to(AP, a, i_s, i_f, g, gto, EK)*np.array(ItoX_dat)/np.array(Ito)/abs(np.mean(ItoX_dat))*abs(np.mean(Ito)))
    ResItoX_t.append(np.array(ItoX_t)-np.array(ItoX_dat)) 
    plt.plot(Ito)
    plt.plot(ItoX_t)
    
    plt.plot(ItoX_dat)
    ItoY_t.append(I_to(np.array(AP), a, i_s, i_f, g, gto, np.array(EK))*np.array(ItoY_dat)/np.array(Ito)/abs(np.mean(ItoY_dat))*abs(np.mean(Ito)))                 
    ResItoY_t.append(np.array(ItoY_t)-np.array(ItoY_dat)) 
#    Isus = I_sus(V, gsus, EK)
    Iks =I_Ks(AP, gks, Ca2p_sl, EK)
    Iks_t.append(I_Ks(AP, gks, Ca2p_sl, EK)*np.array(Iks_dat/Iks/abs(np.mean(Iks_dat))*abs(np.mean(Iks)))                 
    ResIks_t.append(Iks_t-np.array(Iks_dat)) 
    If = I_f(AP, gf_K, gf_Na, EK, ENa)
    If_t.append(I_f(AP, gf_K, gf_Na, EK, ENa)*If_dat/If/abs(np.mean(If_dat))*abs(np.mean(If)))                 
    ResIf_t.append(np.array(If_t)-np.array(If_dat)) 
    Ik1 = I_K1(AP, gk1, CKpout, EK) #inward
    Ik1_t.append(I_K1(AP, gk1, CKpout, EK)*Ik1_dat/Ik1/abs(np.mean(Ik1_dat))*abs(np.mean(Ik1)))                 
    ResIk1_t.append(np.array(Ik1_t)-np.array(Ik1_dat))     
    ICaALL = I_CaL(AP, gCaL)
    ICaALL_t.append(ICaALL[0]*ICaL_dat/ICaALL[0]/abs(np.mean(ICaL_dat))*abs(np.mean(ICaALL[0])))                 
    ResICaALL_t.append(np.array(ICaALL[0])-np.array(ICaL_dat))     
    INa_f = I_Na(AP, h_f, h_s, G_Na, ENa)
    INa_f_t.append(I_Na(AP, h_f, h_s, G_Na, ENa)*INa_dat/INa_f/abs(np.mean(INa_dat))*abs(np.mean(INa_f)))                 
    ResINa_f_t.append(np.array(INa_f_t)-np.array(INa_dat))  
    Ikr = I_Kr(AP, CKpout, xr_s, xr_f, gkr, EK) #outward
    Ikr_t.append(I_Kr(AP, CKpout, xr_s, xr_f, gkr, EK)*Ikr_dat/Ikr/abs(np.mean(Ikr_dat))*abs(np.mean(Ikr)))
    ResIkr_t.append(np.array(Ikr_t)-np.array(Ikr_dat))
    
    I_stim = -40 #muA/muF, use regular, assume amplitude    
    Cm = -(ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+ICaL_dat+INa_dat+Ikr_dat-I_stim+Iothers)/max(dVdt)
    dVdt_1= np.outer(-(Iks_t+ItoX_t+ItoY_t+If_t+Ik1_t+ICaALL_t+INa_f_t+Ikr_t-I_stim+Iothers),1/Cm)
    D = 1.162 #cm**2//=s
    V0_back = Vm
    V =[]
    I = -(ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat+Iothers-I_stim)
    for t in range(201):
        dVdt = FV(np.array(D), V0, np.linspace(0, len(V0), len(V0-1)),I, 0.1)  
        #assume t = 1s 
        dVdt[np.isnan(dVdt)]=0
        V0 +=  dVdt/1000
        print(V0[0])
        V.append(V0)
        I = I_to(V0, a, i_s, i_f, g, gto, EK)*ItoX_dat/ItoX/abs(np.mean(ItoX_dat))*abs(np.mean(Ito))+I_to(V0, a, i_s, i_f, g, gto, EK)*ItoY_dat/Ito/abs(np.mean(ItoY_dat))*abs(np.mean(Ito))+I_Ks(V0, gks, Ca2p_sl, EK)*np.array(Iks_dat)/np.array(Iks)/abs(np.mean(Iks_dat))*abs(np.mean(Iks))+I_f(V0, gf_K, gf_Na, EK, ENa)*If_dat/If/abs(np.mean(If_dat))*abs(np.mean(If))+I_K1(V0, gk1, CKpout, EK)*Ik1_dat/Ik1/abs(np.mean(Ik1_dat))*abs(np.mean(Ik1))+I_CaL(V0, gCaL)[0]*ICaL_dat/ICaALL[0]/abs(np.mean(ICaL_dat))*abs(np.mean(ICaALL[0]))+I_Na(V0, h_f, h_s, G_Na, ENa)*INa_dat/INa_f/abs(np.mean(INa_dat))*abs(np.mean(INa_f))+I_Kr(V0, CKpout, xr_s, xr_f, gkr, EK)*Ikr_dat/Ikr/abs(np.mean(Ikr_dat))*abs(np.mean(Ikr))+Iothers[i]/1000
        ItoX = I_to(V0, a, i_s, i_f, g, gto, EK)*ItoX_dat/ItoX/abs(np.mean(ItoX_dat))*abs(np.mean(ItoX))
        ItoY = I_to(V0, a, i_s, i_f, g, gto, EK)*ItoY_dat/ItoY/abs(np.mean(ItoX_dat))*abs(np.mean(ItoY))        
        I_Ks = I_Ks(V0, gks, Ca2p_sl, EK)*np.array(Iks_dat)/np.array(Iks)/abs(np.mean(Iks_dat))*abs(np.mean(Iks))
        I_f = I_f(V0, gf_K, gf_Na, EK, ENa)*If_dat/If/abs(np.mean(If_dat))*abs(np.mean(If))
        IK1 = I_K1(V0, gk1, CKpout, EK)*Ik1_dat/Ik1/abs(np.mean(Ik1_dat))*abs(np.mean(Ik1))
        I_CaL = I_CaL(V0, gCaL)[0]*ICaL_dat/ICaALL[0]/abs(np.mean(ICaL_dat))*abs(np.mean(ICaALL[0]))
        I_Na = I_Na(V0, h_f, h_s, G_Na, ENa)*INa_dat/INa_f/abs(np.mean(INa_dat))*abs(np.mean(INa_f))
        I_Kr = I_Kr(V0, CKpout, xr_s, xr_f, gkr, EK)*Ikr_dat/Ikr/abs(np.mean(Ikr_dat))*abs(np.mean(Ikr))
        

#    V = -max(AP)*np.exp(np.outer(dVdt,np.linspace(0, 1, 1000)))
#    V = np.array(np.repeat(dVdt/1000, 201)).reshape(np.size(dVdt),201)+np.exp(np.outer(dVdt/1000,np.linspace(0, 1, 201)))
    Varr = np.array(V)
    plt.figure()
    plt.pcolor(Varr)
    plt.title('Vm(1D)_Spatial_Time')
    plt.figure()
    plt.plot(V[0])
    plt.title('t=0, V distributed along 1D')
    plt.figure()
    plt.plot(Varr[:,0])
    plt.title('x=0, Vchange in 201 ms')
    #or        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(125*rv.cdf(x)[23:])    
    plt.plot(100*rv.cdf(x)[23:])  
    plt.plot(124*rv.cdf(x)[23:]) 
    plt.plot(117*rv.cdf(x)[23:])     
    ax.set_xticklabels(('-1','0','1','2','3','4','5'),ha = 'center')
    plt.xlabel('Cycle length(s)')
    plt.ylabel('APD_reduced(ms)')
    V = -max(AP)
  
        plt.figure
        plt.plot(AP)
        plt.plot(AP30)
        plt.plot(AP60)
        plt.plot(AP90)
        plt.legend()
        plt.title('AP reduced in endo, mid and epi')
   
    for t in np.linspace(0,1,999):
        dt = 1/1000
        
        #if len(Vtemp) == 50:
           # spl = sci.interpolate.splrep(np.linspace(0,len(V0),len(V0)), dVdt = FV(D, V0, np.linspace(0, len(V0_X), len(V0_X-1)),ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim, 0.1))
           # dVdt = sci.interpolate.splev(V0[np.linspace(0, t, len(V0))], FV(D, V0[np.linspace(0, t, len(V0))], np.linspace(0, len(V0), len(V0-1)),ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+ICaL_dat+INa_dat+Ikr_dat-I_stim, 0.1))  
        V += dVdt*dt
        Vtemp.append(V)
    plt.plot(Vtemp)
    
    
    
    Ito = []
    ItoX_t = []
    ResItoX_t = []
    ItoY_t = []
    ResItoY_t = []
    I_Ks = []
    Iks_t = []
    ResIks_t = []
    If = []
    If_t = []
    ResIf_t = []
    Ik1 = []
    Ik1_t = []
    ResIk1_t = []
    ICaALL = []
    ICaALL_t = []
    ResICaALL_t = []
    INa_f = []
    INa_f_t = []
    ResINa_f_t = []
    Ikr = []
    Ikr_t = []
    ResIkr_t = []   
    
    AP = np.array(AP)
     EK = np.array(EK)
     ENa = np.array(ENa)
    #activate AP
    Ito = I_to(AP, a, i_s, i_f, g, gto, EK)
    ItoX_t = I_to(AP, a, i_s, i_f, ItoX_dat/Ito*g, ItoX_dat/Ito*gto, EK)
    ItoX_t.append(I_to(AP, a, i_s, i_f, g, gto, EK)*ItoX_dat/Ito/abs(np.mean(ItoX_dat))*abs(np.mean(Ito)))
    ResItoX_t.append(np.array(ItoX_t)-np.array(ItoX_dat)) 
    plt.plot(Ito)
    plt.plot(ItoX_t)
    plt.plot(ItoX_dat)
    ItoY_t.append(I_to(AP, a, i_s, i_f, g, gto, EK)*ItoY_dat/Ito/abs(np.mean(ItoY_dat))*abs(np.mean(Ito)))                 
    ResItoY_t.append(np.array(ItoY_t)-np.array(ItoY_dat)) 
#    Isus = I_sus(V, gsus, EK)
    Iks =I_Ks(np.array(AP), gks, Ca2p_sl, EK)
    Iks_t.append(I_Ks(AP, gks, Ca2p_sl, EK)*np.array(Iks_dat)/np.array(Iks)/abs(np.mean(Iks_dat))*abs(np.mean(Iks)))                 
    ResIks_t.append(np.array(Iks_t)-np.array(Iks_dat)) 
    If = I_f(AP, gf_K, gf_Na, EK, ENa)
    If_t.append(I_f(AP, gf_K, gf_Na, EK, ENa)*If_dat/If/abs(np.mean(If_dat))*abs(np.mean(If)))                 
    ResIf_t.append(np.array(If_t)-np.array(If_dat)) 
    Ik1 = I_K1(AP, gk1, CKpout, EK) #inward
    Ik1_t.append(I_K1(AP, gk1, CKpout, EK)*Ik1_dat/Ik1/abs(np.mean(Ik1_dat))*abs(np.mean(Ik1)))                 
    ResIk1_t.append(np.array(Ik1_t)-np.array(Ik1_dat))     
    ICaALL = I_CaL(AP, gCaL)
    ICaALL_t.append(I_CaL(AP, gCaL)[0]*ICaL_dat/ICaALL[0]/abs(np.mean(ICaL_dat))*abs(np.mean(ICaALL[0])))                 
    ResICaALL_t.append(np.array(Ik1_t)-np.array(Ik1_dat))     
    INa_f = I_Na(AP, h_f, h_s, G_Na, ENa)
    INa_f_t.append(I_Na(AP, h_f, h_s, G_Na, ENa)*INa_dat/INa_f/abs(np.mean(INa_dat))*abs(np.mean(INa_f)))                 
    ResINa_f_t.append(np.array(INa_f_t)-np.array(INa_dat))  
    Ikr = I_Kr(AP, CKpout, xr_s, xr_f, gkr, EK) #outward
    Ikr_t.append(I_Kr(AP, CKpout, xr_s, xr_f, gkr, EK)*Ikr_dat/Ikr/abs(np.mean(Ikr_dat))*abs(np.mean(Ikr)))
    ResIkr_t.append(np.array(Ikr_t)-np.array(Ikr_dat))
    

    
    I_stim = -40 #muA/muF, use regular, assume amplitude    
    Cm = -(ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim)/max(dVdt)
    dVdt_1= np.outer(-(ItoX_t+ItoY_t+Iks_t+If_t+Ik1_t+np.array(ICaALL_t)+INa_f_t+Ikr_t-I_stim),1/Cm)
    D = 1.162 #cm**2//=s
    dVdt = FV(D, IKr0_X, np.linspace(0, len(IKr0_X), len(IKr0_X-1)),ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim, 0.1)  
    #assume t = 1s 
    V = -max(AP)*np.exp(np.outer(dVdt,np.linspace(0, 1, 1000)))
    V_1 = -max(AP)*np.exp(np.outer(dVdt,np.linspace(0, 1, 1000)))
    plt.figure()
    plt.pcolor(V)
    plt.title('ECG(1D)_Spatial_Time')
    V = -max(AP)
    Vtemp = []
    for t in np.linspace(0,1,999):
        dt = 1/1000
        if len(Vtemp) = 50:
            spl = sci.interpolate.splrep(np.linspace(0,len(IKr0_X),len(IKr0_X)), dVdt = FV(D, IKr0_X, np.linspace(0, len(IKr0_X), len(IKr0_X-1)),ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim, 0.1))
            dVdt = sci.interpolate.splev(np.linspace(IKr0_X[np.linspace(0, len(IKr0_X[0:t]), len(IKr0_X))]), FV(D, IKr0_X[0:t], np.linspace(0, len(IKr0_X), len(IKr0_X-1)),ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim, 0.1))  
            V += dVdt*dt
        Vtemp.append(V)
    plt.plot(V)
    
#    XX, YY = np.meshgrid(range(1000), range(201))
 #   plt.figure()
#    plt.contour(XX,YY, V)
#    plt.title()
    V = []
    temp = -max(AP)  
    V.append(temp)
    for t in range(np.size(dVdt)):
        dt = 0.2
        temp += dt*dV[t]
        V.append(temp)
    plt.plot(V)

#    plt.plot(max(AP)*np.exp(np.dot(-dVdt,np.linspace(0,np.size(dVdt)-1))))
    # Bayes Inference on Tpe, tQT, rheart, Tamp, Tarea, TsymG, TsymA
    for t in range(len(AP),len(AP),len(IKr0_X)):
        AP = V[t-len(AP):t]
        V.append(-max(AP))
        ItoX_t.append(I_to(AP, a, i_s, i_f, g, gto, EK)*ItoX_dat/Ito/abs(np.mean(ItoX_dat))*abs(np.mean(Ito)))
        ResItoX_t.append(np.array(ItoX_t)-np.array(ItoX_dat)) 
        ItoY_t.append(I_to(AP, a, i_s, i_f, g, gto, EK)*ItoY_dat/Ito/abs(np.mean(ItoY_dat))*abs(np.mean(Ito)))                 
        ResItoY_t.append(np.array(ItoY_t)-np.array(ItoY_dat)) 
    #    Isus = I_sus(V, gsus, EK)
        Iks_t.append(I_Ks(AP, gks, Ca2p_sl, EK)*Iks_dat/Iks/abs(np.mean(Iks_dat))*abs(np.mean(Iks)))                 
        ResIks_t.append(np.array(Iks_t)-np.array(Iks_dat)) 
        If_t.append(I_f(AP, gf_K, gf_Na, EK, ENa)*If_dat/If/abs(np.mean(If_dat))*abs(np.mean(If)))                 
        ResIks_t.append(np.array(If_t)-np.array(If_dat)) 
        Ik1_t.append(I_K1(AP, gk1, CKpout, EK)*Ik1_dat/Ik1/abs(np.mean(Ik1_dat))*abs(np.mean(Ik1)))                 
        ResIk1_t.append(np.array(Ik1_t)-np.array(Ik1_dat))     
        ICaALL_t.append(I_CaL(AP, gCaL)[0]*ICaL_dat/ICaALL[0]/abs(np.mean(ICaL_dat))*abs(np.mean(ICaALL[0])))                 
        ResICaALL_t.append(np.array(I_CaL[0])-np.array(ICaL_dat))     
        INa_f_t.append(I_Na(AP, h_f, h_s, G_Na, ENa)*INa_dat/INa_f/abs(np.mean(INa_dat))*abs(np.mean(INa_f)))                 
        ResIks_t.append(np.array(INa_f_t)-np.array(INa_dat))  
        Ikr_t.append(I_Kr(AP, CKpout, xr_s, xr_f, gkr, EK)*Ikr_dat/Ikr/abs(np.mean(Ikr_dat))*abs(np.mean(Ikr)))                 
        ResIkr_t.append(np.array(Ikr_t)-np.array(Ikr_dat))
        
    plt.figure()
    plt.plot(ItoX_t)
    plt.title('ItoX_sim')
    plt.figure()
    plt.plot(ResItoX_t)
    plt.title('Res_ItoX')
    plt.figure()
    plt.plot(ItoY_t[:])
    plt.title('ItoY_sim')
    plt.figure()
    plt.plot(ResItoY_t[:])
    plt.title('Res_ItoY')
#    plt.plot(Ito[np.array(np.linspace(0,199*round(len(Ito-1)/200),200), dtype = int)])
    plt.figure()    
    plt.plot(Iks) 
    plt.title('Iks_sim')
    plt.figure()
    plt.plot(ResIks_t)
    plt.title('Res_Iks')    
    plt.figure()
    plt.plot(If_t)
    plt.title('If_sim')
    plt.figure()    
    plt.plot(Ik1_t)
    plt.title('Ik1_sim')   
    plt.figure()
    plt.plot(ICaL[0])
    plt.title('ICaL_sim')     
    plt.figure()    
    plt.plot(INa_f_t) 
    plt.title('INa_sim')      
    plt.figure()    
    plt.plot(Ikr_t) 
    plt.title('Ikr_sim')    
#    Cm = -(ItoX_dat+ItoY_dat+Iks_dat+If_dat+Ik1_dat+np.array(ICaL_dat)+INa_dat+Ikr_dat-I_stim)/dVdt
#    dVdt = -(ItoX_t+ItoY_t+Iks_t+If_t+Ik1_t+np.array(ICaALL_t)+INa_f_t+Ikr_t-I_stim)/Cm
#    dV = dVdt/(max(dVdt)-min(dVdt))
#    plt.figure()
#    plt.plot(dV)
#    plt.title('Activation_Itotal')
    V = []
    temp = -max(AP)  
    V.append(temp)
    for t in range(np.size(dVdt)):
        dt = 0.2
        temp += dt*dV[t]
        V.append(temp)
    plt.plot(V)
    
    xnew = np.linspace(0,len(E1Arr), round(len(E1Arr)*fs/100))
#    xnew = np.linspace(0,len(E1Arr), len(V))
    f = scipy.interpolate.splrep(np.linspace(0,len(E1Arr)-1, len(E1Arr)), E1Arr)
    E1Arr = scipy.interpolate.splev(xnew, f) 
    plt.figure()    
    plt.plot(E1Arr) 
    plt.title('Interpolated_Ikr0')
    
    d = np.sqrt(abs(IKr0_X-IKr0_Y)**2)
    r = 0.0175
    #phi_xy = Extracellular_Unipolar_Potential(CKpin, CKpout,  np.array([IKr0_X, IKr0_Y]), d, r)
#    Vm = np.array([IKr0_X, IKr0_Y])
    Vm = AP
    phi_e = r**2*CKpin/4*CKpout*sum(-np.dot(np.gradient(Vm),np.gradient(1/d)))
    plt.plot(phi_e)
    
    #consider Iks:
    mSample = 10
    Samples = []
    Mus = []
    Sigmas =[]
    p = []
    mu = np.mean(Iks)
    sigma = np.std(Iks)
    #consider ICaL 
    Samples1 = []
    Mus1 = []
    Sigmas1 =[]
    p1 = []
    mu1 = np.mean(ICaL)
    sigma1 = np.std(ICaL)
    for i in range(mSample):
        Samples.append(np.random.normal(mu, sigma, 100))
        Samples1.append(np.random.normal(mu1,sigma1, 100))
    
    Min = min(Iks)
    Max = max(Iks)
 #   mDist = 'PIC50' # or 'HILL'
    mDist = 'HILL'

    
    for i in range(mSample):
    #for Pic 50:
        if mDist == 'PIC50':
            Max = 12
            Min = -12
            Mus.append(Min+i*(Max-Min)/(len(Samples[i])-1))
            sum = 0
            for j in range(mSample):
                sum += Samples[j]-Mus[j]
            Sigmas.append(np.sqrt(sum**2/(len(Samples[i])-1)))
            p.append(np.exp(-(Samples[i]-Mus[i]))/Sigmas[i]/(Sigma[i]*(1+np.exp((Samples[i]-Mus[i])/Sigma[i]))*(1+np.exp(-(Samples[i]-Mus[i])/Sigma[i]))))
        elif mDist == 'HILL':
            Max = 10
            Min = 0.1
            Mus1.append(Min+i*(Max-Min)/(len(Samples1[i])-1))
            sum = 0
            for j in range(mSample):
                sum += Samples1[j]-Mus1[j]
            Sigmas1.append(np.sqrt(sum**2/(len(Samples[i])-1)))
            p1.append(np.exp(-(Samples1[i]-Mus1[i]))/Sigmas1[i]/(Sigma1[i]*(1+np.exp((Samples1[i]-Mus1[i])/Sigma1[i]))*(1+np.exp(-(Samples1[i]-Mus1[i])/Sigma1[i]))))
        else:
            p = np.random()
            Mus.append((Samples[i]+p)/2)
            sum = 0
            for j in range(mSample):
                sum += Samples[j]-Mus[j]
            Sigmas.append(np.sqrt(sum**2/(len(Samples[i])-1)))
            p.append(sum/(len(Samples[i])-1)))

            
            
    #test with one AP_endo_valid:
 
        AP_endo_test = np.real(AP_endo[14000:17100])
        plt.figure(),
        plt.plot(AP_endo_test)
        plt.title('AP_endo_Data')
    
        AP_mid_test = np.real(AP_mid[43800-round(len(E1_cor)/10*5+1):45100-round(len(E1_cor)/10*5+1)])
        plt.figure(),
        plt.plot(AP_mid_test)
        plt.title('AP_mid_Data')  
        AP_epi_test = np.real(AP_epi[1510:2700])
        plt.figure(),
        plt.plot(AP_epi_test)
        plt.title('AP_epi_Data')
    
    h = plt.hist(Samples) 
    plt.title('AP_PIC50_hist')
    h1 = plt.hist(Samples1) 
    plt.title('AP_HILL_hist')    
    Ttest = sci.ttest_ind(Samples, Samples1)
    plt.figure()
    plt.boxplot(Samples)
    plt.title('AP_PIC50_stats')
    plt.figure()
    plt.boxplot(Samples1) 
    plt.title('AP_HILL_hist')    
    plt.figure()    
    plt.boxplot(Ttest)
    plt.title('Ttest of the two different APs')
    
    plt.plot(Samples[0])
    
    
    tp = [10, 9000, 15000, 68600]
    N = 1400
    states = [1, 4, 2, 3]
    drifts = [1100, -4500, -1000, 800]
    #calibrate and validate
    N = 1400
    I0 =  (E1Arr).ravel()
    calib = Calibration(V, fs)
    dVmax = 160
    p1 = 1/(450*dt)
#    p2 = 1/40/dVmax
    p2 = -1/40
    p3 = 1/(500*dt)
    p4 = -1/160
    p5 = 1/(500*dt)
    p6 = 1/120
    p7 = 1/(500*dt)
    p8 = -1/120
    output = Calibrate_Validate(I0, calib, gkr, tp, N, states, drifts, dt, p1, p2, p3, p4, p5, p6, p7, p8, EK)
    IkrO = output[0]
    GkrO = output[1]
    ResiduleO = output[2]
    Ikr = output[0]
    Gkr = output[1]
    Residule = output[2]
    #simulation and retest the parameters for ion channel models
    #activate and recovery (do not consider after using E-4031 to eliminate IKr first doing correction)
    steps = 50
    N1 = 600
    N2 = 1500
    Vsim = []
    Isim = []
    Vsim.append(calib[899:900+N2])   
    Vsim.append(calib[249:250+N1])
    Isim.append(I0[4400:4400+N2])
    Isim.append(I0[100:100+N1])             
    EXP = np.exp(Isim)
    A = np.sort(np.linspace(min(I0), max(I0), N))
    B = np.sort(np.linspace(min(EXP), max(EXP), N)/max(EXP)/10)
#    dVsim = Vsim[1:len(Vsim)]-Vsim[0:len(Vsim)-1]
    dVsim = []
    for i in range(len(Vsim)): 
        dVsim.append(Vsim[i][1:len(Vsim[i])]-Vsim[i][0:len(Vsim[i])-1])                     
    Istim = 100 #pA
#    data is annotated when Cm = -(Ik1+Ito+Isus+I_Kr+I_stim)*dt/(dV+0.00001), I_others = 0
    sc = np.linspace(0.5, 1.5, 5)
    gk1sim = gk1*sc
    gtosim = gto*sc
    gsussim = gsus*sc
    gkrsim = gkr*sc
    gkssim = gks*sc
    
    Ito_0 = I_to(Vsim[0][1:len(Vsim[0])], a, i_s, i_f, g, gto, EK)    
    Isus_0 = I_sus(Vsim[0][1:len(Vsim[0])], gsus, EK)
    Ikr_0 = I_Kr(Vsim[0][1:len(Vsim[0])], CKpout, xr_s, xr_f, gkr, EK) #outward
    Iks_0 =I_Ks(Vsim[0][1:len(Vsim[0])], gks, Ca2p_sl, EK)
    Ik1_0 = I_K1(Vsim[0][1:len(Vsim[0])], gk1, CKpout, EK) #inward
    Ito_1 = I_to(Vsim[1][1:len(Vsim[1])], a, i_s, i_f, g, gto, EK)    
    Isus_1 = I_sus(Vsim[1][1:len(Vsim[1])], gsus, EK)
    Ikr_1 = I_Kr(Vsim[1][1:len(Vsim[1])], CKpout, xr_s, xr_f, gkr, EK) #outward
    Iks_1 =I_Ks(Vsim[1][1:len(Vsim[1])], gks, Ca2p_sl, EK)
    Ik1_1 = I_K1(Vsim[1][1:len(Vsim[1])], gk1, CKpout, EK) #inward    
    Cm0 = -(Ik1_0+Ito_0+Isus_0+np.array(Isim[0]))*dt/(np.array(dVsim[0])+0.00001)
    Cm1 = -(Ik1_1+Ito_1+Isus_1+np.array(Isim[1])+Istim)*dt/(np.array(dVsim[1])+0.00001)
    Lsigma = np.cov([np.ones((np.size(Ik1_0sim))),np.ones((np.size(Ik1_0sim))),np.ones((np.size(Ik1_0sim))),np.ones((np.size(Ik1_0sim)))])
    L = np.random.multivariate_normal([0,0,0,0], Lsigma, (np.size(Ik1_0sim), np.size(Ik1_0sim)))
    Params = []
    idx = 0
    for i in range(len(sc)):
            Ito_0sim = I_to(Vsim[0][1:len(Vsim[0])], a, i_s, i_f, g, gtosim[i], EK)    
            Isus_0sim = I_sus(Vsim[0][1:len(Vsim[0])], gsussim[i], EK)
#            Ikr_0 = I_Kr(Vsim[0][1:len(Vsim[0])], CKpout, xr_s, xr_f, gkrsim[i], EK) #outward
            Iks_0sim =I_Ks(Vsim[0][1:len(Vsim[0])], gkssim[i], Ca2p_sl, EK)
            Ik1_0sim = I_K1(Vsim[0][1:len(Vsim[0])], gk1sim[i], CKpout, EK) #inward
            
#            sci.norm.pdf(np.array(Ito_0sim)) 
            countIto0, centerIto0 = np.histogram(Ito_0sim, bins = 5)
            countIsus0, centerIsus0 = np.histogram(Isus_0sim, bins = 5)
            countIks0, centerIks0 = np.histogram(Iks_0sim, bins = 5)
            countIk10, centerIk10 = np.histogram(Ik1_0sim, bins = 5)
            Mu = [np.mean(Ito_0sim) , np.mean(Isus_0sim), np.mean(Iks_0sim), np.mean(Ik1_0sim)]
            SIGMA = np.cov([Ito_0sim , Isus_0sim, Iks_0sim, Ik1_0sim])
            Distr1 = np.random.normal(np.mean(Ito_0sim), np.std(Ito_0sim), np.size(Ito_0sim))*np.random.normal(np.mean(Isus_0sim), np.std(Isus_0sim), np.size(Isus_0sim))* np.random.normal(np.mean(Iks_0sim), np.std(Iks_0sim), np.size(Iks_0sim))*np.random.normal(np.mean(Ik1_0sim), np.std(Ik1_0sim), np.size(Ik1_0sim))
            Distr2 = np.random.multivariate_normal(Mu, SIGMA, (np.size(Ik1_0sim), np.size(Ik1_0sim)))*np.std(Iks_0sim)*np.std(Ik1_0sim)*np.std(Isus_0sim)*np.std(Ito_0sim)
            Ltemp = np.dot(Distr1,Distr2)
            if np.mean(Ltemp)/np.std(Ltemp) > np.mean(L)/np.std(L):
#            if np.mean(Isim[0] - Ikr_0sim) < np.mean(Res_0sim):
                L = Ltemp 
                idx = i
    Params = [sc[idx], gtosim[idx], gsussim[idx], gkssim[idx], gk1sim[idx]]
    Ikr_0sim = Cm0*(np.array(dVsim[0])+0.00001)/dt-Ito_0sim-Isus_0sim-Iks_0sim-Ik1_0sim
    Cm0sim = Cm0/(Ikr_0sim-Isim[0])
    Ikr_0simcor = Ikr_0sim*Cm0sim/(Ikr_0sim-Isim[0])
    Ikr_0simcorB = Ikr_0simcor-np.mean(Ikr_0simcor)

    plt.figure()
    plt.plot(Isim[0],'r-')
    plt.plot(Ikr_0,'g*')
    plt.plot(Ikr_0simcorB, 'b--')
    plt.legend(['origin', 'initiate', 'first simulation under best parameter'])
    plt.title('gkr:'+str(gkrsim[idx])+' '+'gives the best simulation of Ikr recovery I')
    
    plt.figure()
    plt.plot(abs(Ikr_0 - Isim[0]),'r')
    plt.plot(abs(Ikr_0simcorB - Isim[0])/np.mean(abs(Ikr_0simcorB - Isim[0])), 'g')  
    plt.legend(['Residule of initiate', ' Residule of first simulation under best parameter'])
    plt.title('accordding residule')

    Ito_0simO = Ito_0sim   
    Isus_0simO = Isus_0sim
    Iks_0simO = Iks_0sim
    Ik1_0simO = Ik1_0sim
    Ikr_0simO = Ikr_0sim    

    Res_Ikr_0sim = [] 
    Res_Ito_0sim = []
    Res_Isus_0sim = [] 
    Res_Iks_0sim = []    
    Res_Ik1_0sim = []
    Res_Ikr_0sim.append(np.mean(abs(Ikr_0simO - Isim[0])/np.mean(abs(Ikr_0simO - Isim[0]))))
    Res_Ito_0sim.append(np.mean(abs(Ito_0simO - Ito_0)/np.mean(abs(Ito_0simO - Ito_0))))
    Res_Isus_0sim.append(np.mean(abs(Isus_0simO - Isus_0)/np.mean(abs(Isus_0simO - Isus_0))))
    Res_Iks_0sim.append(np.mean(abs(Iks_0simO - Iks_0)/np.mean(abs(Iks_0simO - Iks_0))))
    Res_Ik1_0sim.append(np.mean(abs(Ik1_0simO - Ik1_0)/np.mean(abs(Ik1_0simO - Ik1_0))))

#    steps = 100
    N = 15
    c = 0
    NIto = np.size(Res_Ito_0sim)
    NIsus = np.size(Res_Isus_0sim)
    NIks = np.size(Res_Iks_0sim) 
    NIkr = np.size(Res_Ikr_0sim) 
    NIk1 = np.size(Res_Ik1_0sim)         
    while (NIk1 <= N | NIto <= N | NIsus <= N | NIks <= N | NIkr <= N) & steps >=0 & c>-20:
        Ik1_0sim_n = Ik1_0sim + np.random.normal(0, 1, np.size(Ik1_0sim))
        Ik1_0sim_sc = Cm0/(Ik1_0sim_n-Ik1_0)
        Ik1_0simcor = Ik1_0sim_n*Ik1_0sim_sc/(Ik1_0sim_n-Iks_0)        
        Ik1_0simcorB = Ik1_0simcor-np.mean(Iks_0simcor)
        if (np.mean(abs(Ik1_0simcorB - Ik1_0)/np.mean(abs(Ik1_0simcorB - Ik1_0))) <= Res_Ik1_0sim[-1])  & (sum(Ik1_0simcorB != Ik1_0sim)>0):
            Res_Ik1_0sim.append(np.mean(abs(Ik1_0simcorB - Ik1_0)/np.mean(abs(Ik1_0simcorB - Ik1_0))))        
            Ik1_0sim = Ik1_0simcorB
            c += 1

        Ito_0sim_n = Ito_0sim + np.random.normal(0, 1, np.size(Ito_0sim))
        Ito_0sim_sc = Cm0/(Ito_0sim_n-Ito_0)
        Ito_0simcor = Ito_0sim_n*Ito_0sim_sc/(Ito_0sim_n-Ito_0)        
        Ito_0simcorB = Ito_0simcor-np.mean(Ito_0simcor)
        if (np.mean(abs(Ito_0simcorB - Ito_0)/np.mean(abs(Ito_0simcorB - Ito_0))) <= Res_Ito_0sim[-1]) & (sum(Ito_0simcorB != Ito_0sim)>0):
            Res_Ito_0sim.append(np.mean(abs(Ito_0simcorB - Ito_0)/np.mean(abs(Ito_0simcorB - Ito_0))))        
            Ito_0sim = Ito_0simcorB
            c  += 1

        Isus_0sim_n = Isus_0sim + np.random.normal(0, 1, np.size(Isus_0sim))
        Isus_0sim_sc = Cm0/(Isus_0sim_n-Isus_0)
        Isus_0simcor = Isus_0sim_n*Isus_0sim_sc/(Isus_0sim_n-Isus_0)        
        Isus_0simcorB = Isus_0simcor-np.mean(Isus_0simcor)
        if (np.mean(abs(Isus_0simcorB - Isus_0)/np.mean(abs(Isus_0simcorB - Isus_0))) <= Res_Isus_0sim[-1]) & (sum(Isus_0simcorB != Isus_0sim)>0):
            Res_Isus_0sim.append(np.mean(abs(Isus_0simcorB - Isus_0)/np.mean(abs(Isus_0simcorB - Isus_0))))        
            Isus_0sim = Isus_0simcorB
            c += 1
        
        Iks_0sim_n = Iks_0sim + np.random.normal(0, 1, np.size(Iks_0sim))
        Iks_0sim_sc = Cm0/(Iks_0sim_n-Iks_0)
        Iks_0simcor = Iks_0sim_n*Iks_0sim_sc/(Iks_0sim_n-Iks_0)        
        Iks_0simcorB = Iks_0simcor-np.mean(Iks_0simcor)
        if (np.mean(abs(Iks_0simcorB - Iks_0)/np.mean(abs(Iks_0simcorB - Iks_0))) <= Res_Iks_0sim[-1])  & (sum(Iks_0simcorB != Iks_0sim)>0):
            Res_Iks_0sim.append(np.mean(abs(Iks_0simcorB - Iks_0)/np.mean(abs(Iks_0simcorB - Iks_0))))        
            Iks_0sim = Iks_0simcorB
            c += 1
            
        Ikr_0sim_n = Ikr_0sim + np.random.normal(0, 1, np.size(Ikr_0sim))
        Ikr_0sim_sc = Cm0/(Ikr_0sim_n-Ikr_0)
        Ikr_0simcor = Ikr_0sim_n*Ikr_0sim_sc/(Ikr_0sim_n-Ikr_0)        
        Ikr_0simcorB = Ikr_0simcor-np.mean(Ikr_0simcor)
        if (np.mean(abs(Ikr_0simcorB - Ikr_0)/np.mean(abs(Ikr_0simcorB - Ikr_0))) <= Res_Ikr_0sim[-1])  & (sum(Ikr_0simcorB != Ikr_0sim)>0):
            Res_Ikr_0sim.append(np.mean(abs(Ikr_0simcorB - Ikr_0)/np.mean(abs(Ikr_0simcorB - Ikr_0))))        
            Ikr_0sim = Ikr_0simcorB
            c += 1
        c -= 1
        NIto = np.size(Res_Ito_0sim)
        NIsus = np.size(Res_Isus_0sim)
        NIks = np.size(Res_Iks_0sim) 
        NIkr = np.size(Res_Ikr_0sim)
       
        steps-= 1
        
        
        

        dVmaxdt = []
        Sensitivity = [] 
 
        Res_ks_0sim.append(abs(Iks_0simcorB - Iks_0)/np.mean(abs(Iks_0simcorB - Iks_0)))        

        Ik1_0sim_n = Ik1_0sim + np.random.normal(0, 1, np.size(Ik1_0sim))
        Ik1_0sim_sc = Cm0/(Ik1_0sim_n-Ik1_0)
        Ik1_0simcor = Ik1_0sim_n*Ik1_0sim_sc/(Ik1_0sim_n-Ik1_0)        
        Ik1_0simcorB = Ik1_0simcor-np.mean(Ik1_0simcor)
        Res_Ik1_0sim.append(abs(Ik1_0simcorB - Ik1_0)/np.mean(abs(Ik1_0simcorB - Ik1_0)))        

        
    Ito_t0 = Ito_0simcorB
    Isus_t0 = Isus_0simcorB
    Iks_t0 = Iks_0simcorB
    Ikr_t0 = Ikr_0simcorB
    Ik1_t0 = Ik1_0simcorB
    dVmax = -(Ik1_t+Ito_t+Isus_t+Ikr_t+Iks_t)/Cm0
    Res_Ik1_t0 = Res_Ik1_0sim[0:14]
    Res_Ito_t0 = Res_Ito_0sim[0:14]
    Res_Isus_t0 = Res_Isus_0sim[0:14]
    Res_Iks_t0 = Res_Iks_0sim[0:14]    
    Res_Ikr_t0 = Res_Ikr_0sim[0:14]
    # Validate 3: DAD
    Ito_DAD = gtosim[idx]*(DAD-EK)    
    Isus_DAD = gsussim[idx]*(DAD-EK)
    Iks_DAD = gkssim[idx]*(DAD-EK)
    Ikr_DAD = gkrsim[idx]*(DAD-EK)
    Ik1_DAD = gk1sim[idx]*(DAD-EK)
    
    Ito_DAD = DAD- Ito_0simO/gtosim[idx]    
    Isus_DAD = gsussim[idx]*(DAD-EK)
    Iks_DAD = gkssim[idx]*(DAD-EK)
    Ikr_DAD = gkrsim[idx]*(DAD-EK)
    Ik1_DAD = gk1sim[idx]*(DAD-EK)
    # Validate 4: EAD
    Ito_EAD = gtosim[idx]*(EAD-EK)
    Isus_EAD = gsussim[idx]*(EAD-EK)
    Iks_EAD = gkssim[idx]*(EAD-EK)
    Ikr_EAD = gkrsim[idx]*(EAD-EK)
    Ik1_EAD = gk1sim[idx]*(EAD-EK)    
    # Validate 4: AP   
    Ito_AP = gtosim[idx]*(AP-EK)    
    Isus_AP = gsussim[idx]*(AP-EK)
    Iks_AP = gkssim[idx]*(AP-EK)
    Ikr_AP = gkrsim[idx]*(AP-EK)
    Ik1_AP = gk1sim[idx]*(AP-EK)

    plt.figure()
    plt.subplot(3,2,1)
    plt.plot(DAD) 
    plt.title('protocal DAD') 
    plt.subplot(3,2,2)
    plt.plot(Ito_DAD)
    plt.ylabel('Gto') 
    plt.subplot(3,2,3)
    plt.plot(Isus_DAD)
    plt.ylabel('Gsus')
    plt.subplot(3,2,4)     
    plt.plot(Iks_DAD)
    plt.ylabel('Gks')
    plt.subplot(3,2,5)
    plt.plot(Ikr_DAD)
    plt.ylabel('Gkr')
    plt.subplot(3,2,6)
    plt.plot(Ik1_DAD)    
    plt.ylabel('Gk1')  
    
    plt.figure()    
    plt.subplot(3,2,1)
    plt.plot(EAD) 
    plt.title('protocal EAD') 
    plt.subplot(3,2,2)
    plt.plot(Ito_EAD)
    plt.ylabel('Gto') 
    plt.subplot(3,2,3)
    plt.plot(Isus_EAD)
    plt.ylabel('Gsus')
    plt.subplot(3,2,4)     
    plt.plot(Iks_EAD)
    plt.ylabel('Gks')
    plt.subplot(3,2,5)
    plt.plot(Ikr_EAD)
    plt.ylabel('Gkr')
    plt.subplot(3,2,6)
    plt.plot(Ik1_EAD)    
    plt.ylabel('Gk1') 

    plt.figure()
    plt.subplot(3,2,1)
    plt.plot(AP) 
    plt.title('protocal AP') 
    plt.subplot(3,2,2)
    plt.plot(Ito_AP)
    plt.ylabel('Gto') 
    plt.subplot(3,2,3)
    plt.plot(Isus_AP)
    plt.ylabel('Gsus')
    plt.subplot(3,2,4)     
    plt.plot(Iks_AP)
    plt.ylabel('Gks')
    plt.subplot(3,2,5)
    plt.plot(Ikr_AP)
    plt.ylabel('Gkr')
    plt.subplot(3,2,6)
    plt.plot(Ik1_AP)    
    plt.ylabel('Gk1') 
    
    
    plt.figure()
    plt.subplot(5,2,1)
    plt.plot(Ito_t0[0:200])
    plt.ylabel('Ito')
    plt.subplot(5,2,2)
    plt.plot(Ito_t0[0:400])
    plt.ylabel('Ito')
    
    
    plt.subplot(5,2,3)
    plt.plot(Isus_t0[0:200])
    plt.ylabel('Isus')
    plt.subplot(5,2,4)
    plt.plot(Isus_t0[0:400])
    plt.ylabel('Isus') 
    
    plt.subplot(5,2,5)
    plt.plot(Iks_t0[0:200])
    plt.ylabel('Iks')     
    plt.subplot(5,2,6)
    plt.plot(Iks_t0[0:400])
    plt.ylabel('Iks')   
    
    plt.subplot(5,2,7)
    plt.plot(Ikr_t0[0:200])
    plt.ylabel('Ikr')     
    plt.subplot(5,2,8)
    plt.plot(Ikr_t0[0:400])
    plt.ylabel('Ikr')
    
    plt.subplot(5,2,9)
    plt.plot(Ik1_t0[0:200])
    plt.ylabel('Ik1')     
    plt.subplot(5,2,10)
    plt.plot(Ik1_t0[0:400])
    plt.ylabel('Ik1')
    
    
    plt.figure()
    plt.subplot(5,3,1)
    plt.plot(Ito_t0[0:200])
    plt.ylabel('Ito')
    plt.subplot(5,3,2)
    plt.plot(Ito_t0[0:200])
    plt.ylabel('Ito')     
    
    plt.subplot(5,3,3)
    plt.plot(Isus_t0[0:200])
    plt.ylabel('Isus')     
    plt.subplot(5,3,4)
    plt.plot(Isus_t0[0:400])
    plt.ylabel('Isus') 
    
    plt.subplot(5,2,5)
    plt.plot(Iks_t0[0:200])
    plt.ylabel('Iks')     
    plt.subplot(5,2,6)
    plt.plot(Iks_t0[0:400])
    plt.ylabel('Iks')   
    
    plt.subplot(5,2,7)
    plt.plot(Ikr_t0[0:200])
    plt.ylabel('Ikr')     
    plt.subplot(5,2,8)
    plt.plot(Ikr_t0[0:400])
    plt.ylabel('Ikr')
    
    plt.subplot(5,2,9)
    plt.plot(Ik1_t0[0:200])
    plt.ylabel('Ik1')     
    plt.subplot(5,2,10)
    plt.plot(Ik1_t0[0:400])
    plt.ylabel('Ik1')    
    
    

#    dVmaxdt =       

    mean_Res_Ik1 = []
    mean_Res_Ito = []
    mean_Res_Isus = []  
    mean_Res_Iks = []
    mean_Res_Ikr = []
    std_Res_Ik1 = []
    std_Res_Ito = []
    std_Res_Isus = []  
    std_Res_Iks = []
    std_Res_Ikr = []
    CV_Res_Ik1 = []
    CV_Res_Ito = []
    CV_Res_Isus = []  
    CV_Res_Iks = []
    CV_Res_Ikr = []
    
    for i in range(1,15):
        mean_Res_Ik1.append(np.mean(Res_Ik1_t0[0:i]))
        mean_Res_Ito.append(np.mean(Res_Ito_t0[0:i]))
        mean_Res_Isus.append(np.mean(Res_Isus_t0[0:i]))  
        mean_Res_Iks.append(np.mean(Res_Iks_t0[0:i]))
        mean_Res_Ikr.append(np.mean(Res_Ikr_t0[0:i]))
        std_Res_Ik1.append(np.std(Res_Ik1_t0[0:i]))
        std_Res_Ito.append(np.std(Res_Ito_t0[0:i]))
        std_Res_Isus.append(np.std(Res_Isus_t0[0:i])) 
        std_Res_Iks.append(np.std(Res_Iks_t0[0:i]))
        std_Res_Ikr.append(np.std(Res_Ikr_t0[0:i]))
        CV_Res_Ik1.append(mean_Res_Ik1[-1]/(std_Res_Ik1[-1]+0.000001))
        CV_Res_Ito.append(mean_Res_Ito[-1]/(std_Res_Ito[-1]+0.000001))
        CV_Res_Isus.append(mean_Res_Isus[-1]/(std_Res_Isus[-1]+0.000001)) 
        CV_Res_Iks.append(mean_Res_Iks[-1]/(std_Res_Iks[-1]+0.000001))
        CV_Res_Ikr.append(mean_Res_Ikr[-1]/(std_Res_Ikr[-1]+0.000001))        
        

    Sensitivity = Sensitivity(np.array([Ito_0, Isus_0, Iks_0, Ik1_0, Ikr_0]), np.array([Ito_0simO, Isus_0simO, Iks_0simO, Ik1_0simO, Ikr_0simO]))


    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(mean_Res_Ik1) 
    plt.title('Ik1')
    plt.ylabel('Res_mean')
    plt.subplot(2,2,2)
    plt.plot(std_Res_Ik1)
    plt.ylabel('Res_std') 
    plt.subplot(2,2,3)
    plt.plot(CV_Res_Ik1)
    plt.ylabel('Res_CV')
    plt.subplot(2,2,4)     
    plt.plot(Sensitivity[3,:])
    plt.ylabel('sensitivity')
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(mean_Res_Ito) 
    plt.title('Ito')
    plt.ylabel('Res_mean')
    plt.subplot(2,2,2)
    plt.plot(std_Res_Ito)
    plt.ylabel('Res_std') 
    plt.subplot(2,2,3)
    plt.plot(CV_Res_Ito)
    plt.ylabel('Res_CV')
    plt.subplot(2,2,4)     
    plt.plot(Sensitivity[0,:])
    plt.ylabel('sensitivity')

    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(mean_Res_Iks) 
    plt.title('Iks')
    plt.ylabel('Res_mean')
    plt.subplot(2,2,2)
    plt.plot(std_Res_Iks)
    plt.ylabel('Res_std') 
    plt.subplot(2,2,3)
    plt.plot(CV_Res_Iks)
    plt.ylabel('Res_CV')
    plt.subplot(2,2,4)     
    plt.plot(Sensitivity[2,:])
    plt.ylabel('sensitivity')
    
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(mean_Res_Ikr) 
    plt.title('Isus')
    plt.ylabel('Res_mean')
    plt.subplot(2,2,2)
    plt.plot(std_Res_Isus)
    plt.ylabel('Res_std') 
    plt.subplot(2,2,3)
    plt.plot(CV_Res_Isus)
    plt.ylabel('Res_CV')
    plt.subplot(2,2,4)     
    plt.plot(Sensitivity[1,:])
    plt.ylabel('sensitivity')

    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(mean_Res_Ikr) 
    plt.title('Ikr')
    plt.ylabel('Res_mean')
    plt.subplot(2,2,2)
    plt.plot(std_Res_Ikr)
    plt.ylabel('Res_std') 
    plt.subplot(2,2,3)
    plt.plot(CV_Res_Ikr)
    plt.ylabel('Res_CV')
    plt.subplot(2,2,4)     
    plt.plot(Sensitivity[4,:])
    plt.ylabel('sensitivity')

    muIkr_0sim = np.mean(Ikr_0sim)
    sigmaIkr_sim = np.std(Ikr_0sim)
#    plt.plot(Res_0sim,sci.norm.interval(alpha=0.05) )
    count0, center0 = np.histogram(Ikr_0sim)
    p0 = 1/1+np.exp(-count0/sum(count0))

    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(Ik1_0sim, density = True)
    plt.title('Ik1')
    plt.subplot(2,1,2)
    plt.hist(Sensitivity[3,], density = True)
    plt.title('Ik1 Sensitivity')
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(Ito_0sim, density = True)
    plt.title('Ito')
    plt.subplot(2,1,2)
    plt.hist(Sensitivity[0,], density = True)
    plt.title('Ito Sensitivity')  
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(Iks_0sim, density = True)
    plt.title('Iks')
    plt.subplot(2,1,2)
    plt.hist(Sensitivity[2,], density = True)
    plt.title('Iks Sensitivity')
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(Ikr_0sim, density = True)
    plt.title('Ikr')
    plt.subplot(2,1,2)
    plt.hist(Sensitivity[4,], density = True)
    plt.title('Ikr Sensitivity')   
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.hist(Isus_0sim, density = True)
    plt.title('Isus')
    plt.subplot(2,1,2)
    plt.hist(Sensitivity[1,], density = True)
    plt.title('Isus Sensitivity')      
    
    
    #EAD, BAD
    
    
    
    
    