# =============================================================================
# # -*- coding: utf-8 -*-
#Author:Lilyheasamiko
# =============================================================================
"""
Spyder Editor

This is a temporary script file.
"""

#CML_noble_varghese_kohl_noble_1998_basic_with_sac
import numpy as np
import matplotlib.pyplot as plt
#O stands for dopamine hydroxide
#R stands for reduced after phosphate buffer solution (PB)
def PeakCurrentIntensity(C, alpha, n):
    A = 0.141#electrode area(cm**2)
    D = 0.67*10**(-5) #diffusion coefficient of the analyte(cm**2/s)
    mu0 = 1#scan rate(V/s) 0.93
    ip =2.99*10**5*n**3/2*alpha**0.5*A*C*D**0.5*mu0**0.5
    return ip

def PeakCurrentIntensity(C, alpha, n):
    A = 0.141#electrode area(cm**2)
    D = 0.67*10**(-5) #diffusion coefficient of the analyte(cm**2/s)
    mu0 = 1#scan rate(V/s) 0.93
    ip =2.99*10**5*n**3/2*alpha**0.5*A*C*D**0.5*mu0**0.5
    return ip

def roughness(A):
    Ag = 0.126 # geometry SPCE
    f = A/Ag # roughness factor  
    return f

def SPEC(alpha, u, n):
    k = 9.9*10**(-3)# heterogeneous rate constant(cm/s)
    D0 = 0.67*10**(-5)# diffusion coefficient for the dopamine
    DR = 1*10**(-4)# diffusion coefficient for dopaminequinone (approximated to the dopamine diffusion coefficient) (cm2/s)
    R = 8.314 # unniversal gas constant(J/mol K)
    T = 273+27.15 # room temperature
    F = 96485.3329 #Faraday constant
    psi= k*(D0/DR)**(alpha/2)*(R*T)**0.5*(np.pi*n*F*D0*u)**(-0.5)
    return psi

def CapacitiveCurrent(A, u):
    Cdl = 85 #capacitance (F/cm**2)    
    ic = A*Cdl*u
    return ic

def dE(psi):
    dE = np.exp((3.69 - np.log(psi))/1.16)+59
    return dE

def main():
    C = 125 #bulk concentration (mol/cm**3)
    mu = [10, 25, 50, 100]
    print("scan rate:")
    print(mu)
    A = 0.141#electrode area(cm**2)  
    f = roughness(A)
    print("f:")
    print(f)
    k = 9.9*10**(-3)# heterogeneous rate constant(cm/s)
    D0 = 0.67*10**(-5)# diffusion coefficient for the dopamine
    DR = 1*10**(-4)# diffusion coefficient for dopaminequinone (approximated to the dopamine diffusion coefficient) (cm2/s)
    R = 8.314 # unniversal gas constant(J/mol K)
    T = 273+27.15 # room temperature
    F = 96485.3329 #Faraday constant
    alpha = 0.5
    # variables in time
    PSI = []
    IPS = []
    sensitivity = []
    residule =[]
    tCritic = []
    TCritic = []
    Slope = []
    for u in mu:        
        X= []
        for t in range(1, 1000):
            i0 = CapacitiveCurrent(A, u)/100
            Tp = f*u*t 
            CO = C/(1+np.exp(-i0))
            CR = C/(1+np.exp(i0))
            uO = CO/C
            uR = CR/C
            trA = k/(f*u*D0)**0.5
            n = f*u*t#number of electron transferred
            ips = PeakCurrentIntensity(C, alpha, n)/100*10**(-9)
            d = DR/D0 
            x = d*(f*u/D0)*100*10**(-9)
            X.append(ips/i0*x)
        tcritic = np.argmin(abs(np.log(np.array(X))))
        tCritic.append(tcritic)
        residule.append(min(abs(np.log(np.array(X)))))
        ncritic = f*u*tcritic#number of electron transferred
        Tpcritic = f*u*tcritic 
        TCritic.append(Tpcritic)
        sensitivity.append(uO**2+d**0.5*uR**2-(1+d**0.5*np.exp(-PeakCurrentIntensity(C, alpha, ncritic)/100*10**(-9))/(1+np.exp(-PeakCurrentIntensity(C, alpha, ncritic)/100*10**(-9)))))
        Slope.append(np.log(np.exp(f*np.log(u)))/np.log(u))
        Ips = []
        Psi = []
        nt0 = f*u*np.array(range(1, int(round(Tpcritic))))
#        t0 = np.argmin(abs(trA - (PeakCurrentIntensity(C, alpha, nt0)/100*10**(-9))/i0))
        t0 = np.argmin(SPEC(alpha, u, nt0)*100/(trA - (PeakCurrentIntensity(C, alpha, nt0)/100*10**(-9))/i0))
#        dPsi=SPEC(alpha, u, nt0[t0])*100*trA-SPEC(alpha, u, nt0[t0-1])*100
#        dIps = PeakCurrentIntensity(C, alpha, nt0[t0])**0.5*10**(-3)*(trA-1)/nt0[t0]
#        for t in range(1, int(round(Tpcritic))-1):
#            if t <= int(round(Tpcritic))-1-t0:
        for t in range(1, 2*int(round(Tpcritic))-2*t0-2):
            if t < int(round(Tpcritic))-t0-1:
                nt = nt0[int(round(Tpcritic))-1-t] #number of electron transferred
#                Ips.append(-1.1090+np.log(DR/D0*t/100*10**(-9))**0.5)
                Ips.append(PeakCurrentIntensity(C, alpha, f*u*t)**0.5*10**(-3)/(f*u*t)/i0)
                Psi.append(-SPEC(alpha, u, nt)*100)
        #    elif t == t0:
        #        nt = nt0[t0]
#                Ips.append(-0.20802+np.log(trA*t/100*10**(-9)/(np.pi*alpha)**0.5))
        #        Psi.append(SPEC(alpha, u, nt)*100*trA-dPsi)                
        #        Ips.append(PeakCurrentIntensity(C, alpha, nt0[t])**0.5*10**(-3)/nt0[t]/i0)
            else:
        #        nt = nt0[t-146]
                nt = nt0[t-(int(round(Tpcritic))-t0-1-t0)]
#                print(int(round(Tpcritic))-t0-1-t0)
                if t == int(round(Tpcritic))-t0-1:
#                   Ips.append(-0.20802+np.log(trA*t/100*10**(-9)/(np.pi*alpha)**0.5))
                    dPsi = -SPEC(alpha, u, nt)*100*2*trA-Psi[-1]
#                    ntemp =nt0[int(round(Tpcritic))-1-t+1]
#                    dPsi = SPEC(alpha, u, nt)*100*trA-SPEC(alpha, u, ntemp)*100
                    dIps = PeakCurrentIntensity(C, alpha, f*u*t)**0.5*10**(-3)/(f*u*t)/i0-Ips[-4]
                Psi.append(-SPEC(alpha, u, nt)*100*2*trA-dPsi)              
                Ips.append(Ips[-1]+dIps)
        IPS.append(Ips)
        PSI.append(Psi)
        

    print("Slope")
    print(Slope[0])
    
    MeantCritic = []
    MeantCritic.append(np.mean(tCritic))
    MeantTCritic = []
    MeantTCritic.append(np.mean(TCritic))   
    Meanresidule = []
    Meanresidule.append(np.mean(residule))
    Meansensitivity = []
    Meansensitivity.append(np.mean(sensitivity))
    StdtCritic = []
    StdtCritic.append(np.std(tCritic))
    StdTCritic = []
    StdTCritic.append(np.std(TCritic))
    Stdresidule = []
    Stdresidule.append(np.std(residule))
    Stdsensitivity = []
    Stdsensitivity.append(np.std(sensitivity))    
   

    IPSsc = []
    PsiU = []
    IpsU = []
    plt.figure()     
    for i in range(len(mu)):
        sc = max(np.array(IPS)[0])/max(np.array(IPS)[i])
        xx = sc*np.array(IPS[i])
        PsiS = -0.02+0.01*(len(mu)-i)
        IpsS = 4
        IPSsc.append(xx)
        PsiU.append(2*PsiS-np.array(PSI[i]))
        IpsU.append(2*IpsS-xx) 
        plt.plot(IPSsc[i], PSI[i], IpsU[i], PsiU[i])

        if i <  len(mu):
            MeantCritic.append(np.mean(tCritic[0:i+1]))
            MeantTCritic.append(np.mean(TCritic[0:i+1]))   
            Meanresidule.append(np.mean(residule[0:i+1]))
            Meansensitivity.append(np.mean(sensitivity[0:i+1]))
            StdtCritic.append(np.std(tCritic[0:i+1]))
            StdTCritic.append(np.std(TCritic[0:i+1]))
            Stdresidule.append(np.std(residule[0:i+1]))
            Stdsensitivity.append(np.std(sensitivity[0:i+1]))  
        

    plt.figure()
    x = np.array(range(len(MeantCritic)-1))
    y = MeantCritic[1:len(MeantCritic)]
    yerr = StdtCritic[1:len(MeantCritic)]
    plt.errorbar(x, y, yerr[::-1])
    plt.title('criticle transfer time')
    
    plt.figure()
    y = MeantTCritic[1:len(MeantCritic)]
    yerr = StdTCritic[1:len(MeantCritic)]
    plt.errorbar(x, y, yerr[::-1])
    plt.title('period')    
    
    plt.figure()
    y = Meanresidule[1:len(MeantCritic)]
    yerr = Stdresidule[1:len(MeantCritic)]
    plt.errorbar(x, y, yerr[::-1])    
    plt.title('residule')    
    
    plt.figure()
    y = Meansensitivity[1:len(MeantCritic)]
    yerr = Stdsensitivity[1:len(MeantCritic)]
    plt.errorbar(x, y, yerr[::-1])    
    plt.title('sensitivity of the residule')


  
    print("Slope")
    print(Slope)
    print("-----------------------------------------------")
    print(np.log(mu))
    


        