# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 01:11:24 2020

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


class Pts:
    def __init__(self,*args,**kwargs):
        self.BT = kwargs['BT']
        self.DT = kwargs['DT']
        self.IDX = kwargs['IDX']
        self.name = kwargs['name']
        
    def w_copy(self):
        return Pts(BT =self.BT,
                   DT = self.DT,
                   IDX = self.IDX,
                   name = self.name)

    
def CrossingCount(signal, idx0, BT, DT, nQRS, nNoise, PK, PKID):
    count = 0
    pts = []
    for idx in range(len(signal)-1):
        if (signal[idx] - DT)*(DT - signal[idx+1]) >= 0 + (signal[idx] - BT)*(BT - signal[idx+1]) >=0:            
            count += 1
            if count >4:
                pts.append(Pts(BT = 1.5* BT,
                               DT = max(0.5*max(signal), 1.5* BT),
                               IDX = [],
                               name = 'Noise'))
#                return Pts
                print(Pts)
            elif count >= 2 & count <= 4 & idx == len(signal)-1:
                IDX = []
                IDX.append(idx0)   
                IDX.append(PKID)
                IDX.append(idx+idx0+1)
                pts.append(Pts(BT = 0.75*BT + 0.2*max(signal),
                               DT = max(0.5*max(signal), 0.75*BT + 0.2*max(signal)),
                               IDX =IDX, 
                               name = 'QRS'))
#                return Pts
                print(Pts)
            elif idx == len(signal)-1 & count < 2:
                IDX = []
                IDX.append(idx0)    
                IDX.append(idx+idx0+1) 
                pts.append(Pts(BT = 0.5* BT,  
                               DT = 0.5* BT, 
                               IDX = IDX,
                               name = 'DK'))
                print(Pts)
#                return Pts
        elif idx == len(signal)-1:
            if count < 2:
                IDX = []
                IDX.append(idx0)    
                IDX.append(idx+idx0+1) 
                pts.append(Pts(BT = 0.5* BT,           
                               DT = 0.5* BT,   
                               IDX = IDX,
                               name = 'DK'))
                print(Pts)
#                return Pts
            elif count >= 2 & count <= 4:
                IDX = []
                IDX.append(idx0)   
                IDX.append(PKID)
                IDX.append(idx+idx0+1)
                pts.append(Pts(BT = BT,
                               DT = max(0.75*DT, 0.5*BT),
                               IDX = IDX,
                               name = 'QRS'))
                print(Pts)
#                return Pts

def EnergyCrossingCount(signal, Pts, nQRS, nNoise, PK, PKID):    
    lag = 2500/20*96/60 #around 96ms avoid Twave
    lag1 = 2500/20
    ET = []
    for idx in range(len(signal)):
        if idx <= lag1:
            ET.append(4*(0.75*signal[idx]+0.5*PK))
        elif ET <= lag:
            ET.append(0.5*signal[idx])
        else:
            ET.append(0.2*signal[idx])
    for idx in range(len(signal)):
        if sum(np.array(ET[idx:]) - signal[idx:] <0)==0:
            Pts.IDX.append(idx)
            Pts.name = 'QRS'
            return Pts
    return False

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
    P_array = np.array(P_a)
    Q_array = np.array(Q_a)    
#    P_param = np.mean(P_array[random.sample(range(100),20),:],0)  # q(theta|X), p(theta)  
    P_param =  np.mean(P_array,0)
    P_measure = P_array
    KLtemp = P_param*np.log(P_param-P_measure)
    KLtemp[np.isnan(KLtemp)]=0
    KL = sum(KLtemp, 0)
    E = np.mean(P_measure,0)# E[P(X|theta)]
    Obj = -KL+E
    b = np.argmin(Obj)
    Ib = np.linspace(1.1, 2.1, 10)
    Ib_opt = Ib[b]
    P = P_array[:, b]
    Q = Q_array[:, b]    
    return  [P, Q, Ib_opt, Obj] 
    """"
    plt.figure()
    plt.plot(-KL)
    plt.title('-KL')
    plt.figure()   
    plt.plot(E)
    plt.title('Expectation of measurement')
    plt.figure()     
    plt.plot(Obj)
    plt.title('Optimization Object Function')
    """"


def draw2Dand3D(P, Q, source0, dsource, Ib_opt, V, Is, d, t, z, theta, dt, n_ECG):
    source0 = [3.33,6.67,6.67, 0, 0, 0] 
    dsource = [0.3, 1.9, 2.5, -3.7, 3.2, -34.7]
    #mesh 1.67*1.67*4, from source0 = (3.33,6.67,6.67, 0, 0, 0) to (60, 60, 60, ), dsource=(0.3, 1.9, 2.5, -3.7, 3.2, -34.7)
    #t = 10
    #z = z0
    #theta = 
    #iters = 100
    #dz = 0
    dt = 0.5
    dV = V[t:int(round(t+n_ECG/20*dt))]

    zl = 4
    xl = 1.67
    yl = 1.67
    z0 = source0[2]
    x0 = source0[0]
    y0 = source0[1]
    alpha0 = source0[3]
    beta0 = source0[4] 
    theta0 = source0[5]
    dx = dsource[3]/100*x0
    dy = dsource[4]/100*y0
    dz = dsource[5]/100*z0 
    dalpha = dsource[0]/180*np.pi
    dbeta = dsource[1]/180*np.pi
    dtheta = dsource[2]/180*np.pi    
    
    #at fixed time
    #assume angle unchanged
    nx = int(round(xl/abs(dx)+0.5))
    ny = int(round(yl/abs(dy)+0.5))
    Vtxy = np.repeat(0.1, (ny-1)*(nx-1)).reshape(ny-1, nx-1)
    V0 = np.mean(Q)
    for i in range(ny-1):
        for j in range(nx-1):
            Vtxy[i, j] = V0*dx*(i+1)*dy*(j+1)
    plt.figure()
    plt.pcolor(Vtxy)
    #include angle change unit time dt = 0.5s
    dt = 0.5
#    dV = V[t:int(round(t+n_ECG/20*dt))]
    dxs = []
    dys = []
    for i in range(len(dV)-1):
        dxs.append(dV[i+1]/dV[i]*dsource[3])
        dys.append(dV[i+1]/dV[i]*dsource[4])
        dx = np.mean(dxs)
        dy = np.mean(dys)
        STDx = np.std(dxs)
        STDy = np.std(dys)    
        nx1 = int(round(xl/abs(dx*np.cos(alpha0+dalpha)/np.cos(alpha0))+0.5))
        ny1 = int(round(yl/abs(dy*np.cos(beta0+dbeta)/np.cos(beta0))+0.5))
        Vtxy1 = np.repeat(0.1, ny1*nx1).reshape(ny1, nx1)
        V01 = np.mean(Q)
        if ny1>1 and nx1 >1:
            for i in range(ny1-1):
                for j in range(nx1-1):
                    Vtxy1[i, j] = V01*dx*(i+1)*dy*(j+1)
        else:
            #consider neuman neighbor of the pixel(4 neibors)
            nx1 = int(round(xl*3/abs(dx*np.cos(alpha0+dalpha)/np.cos(alpha0))+0.5))
            ny1 = int(round(yl*3/abs(dy*np.cos(beta0+dbeta)/np.cos(beta0))+0.5))

            Vtxy1 = np.repeat(0.1, (ny1*2-1)*(nx1*2-1)).reshape((ny1*2-1), (nx1*2-1))
            for i in range(ny1*2-1):
                for j in range(nx1*2-1):
                    if j < 1 & i == 1:
                        Vtxy1[i, j] = V01*dx*np.cos(alpha0+dalpha)/np.cos(alpha0)*(j+1) 
                    if j == 1 & i <= 1:   
                        Vtxy1[i, j] = V01*dy*np.cos(beta0+dbeta)/np.cos(beta0)*(ny1-i) 
                    if (j == 1 & i == 2) | (j > 1 & i == 1):  
                        Vtxy1[i, j] = V01*dx*np.cos(alpha0+dalpha)/np.cos(alpha0)*(i+1)*dy*np.cos(beta0+dbeta)/np.cos(beta0)*(j+1) 
                    else:
                        Vtxy1[i, j] = 0
            x, y = np.mgrid[range(ny1*2-1), range(nx1*2-1)]
            x_new, y_new = np.mgrid[range(ny1*4), range(nx1*4)]
            scipy.interpolate.griddata((x,y), Vtxy1, (x_new, y_new), method = 'cubic')
                        
    plt.figure()
    plt.pcolor(Vtxy1)

    #deeper under fat including angle change
    #dz = 0
    #dt = 0.5
    xl = 16.67
    
    dalpha = 5.51/180*np.pi
    dx = 13.32/100*x0
    dy = dx/x0*y0
    
    nx1 = int(round(xl/abs(dx*np.cos(alpha0+dalpha)/np.cos(alpha0))+0.5))
#    ny1 = nx1

    Vtxy1 = np.repeat(0.1, (nx1-1)*(nx1-1)).reshape(nx1-1, nx1-1)
    V01 = np.mean(Q)
    for i in range(nx1-1):
        for j in range(nx1-1):
            Vtxy1[i, j] = V01*dx*(i+1)*dy*(j+1)
    plt.figure()
    plt.pcolor(Vtxy1)
    
    #deeper under fat including muscle change
    xl1 = 33.34
    
    dalpha = 12.57/180*np.pi
    dx = 29.45/100*x0
    dy = dx/x0*y0

    nx11 = int(round(xl1/abs(dx*np.cos(alpha0+dalpha)/np.cos(alpha0))+0.5))
    Vtxy11 = np.repeat(0.1, (nx11-1)*(nx11-1)).reshape(nx11-1, nx11-1)
    V01 = np.mean(Q)
    for i in range(nx11-1):
        for j in range(nx11-1):
            Vtxy11[i, j] = V01*dx*(i+1)*dy*(j+1)
    plt.figure()
    plt.pcolor(Vtxy11)

    #consider the artrial
    #for the inside circle on the plane z = z0: x**2+y**2=3.33**2+6.67**2
    zl = 4
    xl = 1.67
    yl = 1.67
    z0 = source0[2]
    x0 = source0[0]
    y0 = source0[1]
    alpha0 = source0[3]
    beta0 = source0[4] 
    theta0 = source0[5]
    dx = dsource[3]/100*x0
    dy = dsource[4]/100*y0
    dz = dsource[5]/100*z0 
    dalpha = dsource[0]/180*np.pi
    dbeta = dsource[1]/180*np.pi
    dtheta = dsource[2]/180*np.pi    
    
    #at fixed time
    #assume angle unchanged

    R2 = 3.33**2+6.67**2
    R = R2**0.5
    cy = []
    for cx in np.linspace(-R, R, int(round((R*2-dx)/abs(dx)+0.5))):
        cy.append((R2-cx**2)**0.5)
    cx =  np.linspace(-R, R, int(round((R*2+dx)/abs(dx)+0.5)))
    cy = np.array(cy)
    idxy = np.meshgrid(cx, cy)
    nry = int(round((max(cy)*2/abs(dy)+1)))
    nrx = int(round((max(cx)*2/abs(dx)+0.5)))
    nrxl = int(round((3.33+R)/abs(dx)))
    nrxr = nrx-nrxl
    nryd = int(round((3.33+R)/abs(dy)))
    nryu = nry-nryd
    Vct = np.repeat(0.1, nrx*nry).reshape(nry, nrx)
    for i in range(nry):
        if i <= nryd:
            for j in range(nrxl):
                Vct[i, j] = V0*abs(dx)*(nryd-i)*dy*(nrxl-j)
            for j in range(nrxl, nrx):
                Vct[i, j] = V0*abs(dx)*(nryd-i)*dy*(j-nrxl)
        else:
            for j in range(nrxl):
                Vct[i, j] = V0*abs(dx)*(i-nryd)*dy*(nrxl-j)
            for j in range(nrxl, nrx):
                Vct[i, j] = V0*abs(dx)*(i-nryd)*dy*(j-nrxl)
    xnew, ynew = np.mgrid[range(-int(round(nry/2)-0.5), int(round(nry/2)+0.5)), range(-int(round(nry/2)-0.5), int(round(nry/2)+0.5))]
    x, y = np.mgrid[range(-int(round(nry/2)-0.5), int(round(nry/2)+0.5)), range(-int(round(nrx/2)-0.5), int(round(nrx/2+0.5)+0.5))]
#    tck = scipy.interpolate.griddata(np.meshgrid(range(nry), range(nrx)), Vct, (xnew, ynew), method ='linear')
    tck = scipy.interpolate.bisplrep(x, y, Vct, s =0)#
    VctNew = scipy.interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    plt.figure()
    plt.pcolor(VctNew)
    plt.plot((cx+8)/2,(R2/4-cx**2/4)**0.5+4, (cx+8)/2,-(R2/4-cx**2/4)**0.5+4)
    plt.plot((3.33+8)/2,(R2/4-3.33**2/4)**0.5+4,'ro')
    plt.title('2D')
    
    #3D consider:z = y**2-6.67**2: x**2+y**2 =3.33**2+(6.67-z)**2 
    cz = np.array(cy**2)-6.67**2
    nz =  int(round(max(cz)/abs(dz)+0.5))
    nzNew =  11
#    cznew = cz[np.array(np.linspace(0,120,11), int)]
    cznew = []
    for i in np.array(np.linspace(0,120,11), int):
         cznew.append(cz[i])
#    cznew = cz[0:]
#    R2 = 3.33**2+(6.67**2+z)**2 
#    R = R2**0.5

    cy = []
    cx = []
    nry = []
    nrx = []
    nrxl = []
    nrxr = []
    nryd = []
    nryu = []
    R1 = []
    
    for cz in cznew:
        R2 = 3.33**2+(6.67**2+cz)**2 
        R = R2**0.5
        R1.append(R)
        tempy = []
        tempx = []
        for x in np.linspace(-R, R, int(round((R*2-dx)/abs(dx)+0.5))):
            tempy.append(abs(R2-x**2)**0.5)
            tempx.append(x)
  #      xnew = np.linspace(-int(round(len(tempx)/2)-0.5), int(round(len(tempx)/2)+0.5), len(tempx)*20)
        xnew = np.linspace(min(tempx), max(tempx), len(tempx)*20)
        ynew = np.linspace(-int(round(len(tempy)/2)-0.5), int(round(len(tempx)/2)+0.5), len(tempy)*20)
        xx = xnew
        f = scipy.interpolate.interp1d(tempx, (abs(R2-np.array(tempx, dtype = float)**2)**0.5))
        yy = f(xnew)
#        yy = scipy.interpolate.griddata(x, ,  , method ='linear')

        cy.append(np.array(yy, dtype = float))
        cx.append(np.array(xx, dtype = float))
        nry.append(int(round((max(tempy)*2/abs(dy)+1))))
        nrx.append(int(round((max(tempx)*2/abs(dx)+0.5))))
        nrxl.append(int(round((3.33+R)/abs(dx))))
        nrxr.append(nrx[-1]-nrxl[-1])
        nryd.append(int(round((3.33+R)/abs(dy))))
        nryu.append(nry[-1]-nryd[-1])
    cy = np.array(cy)
    cx = np.array(cx)
    idxy = np.meshgrid(cx, cy)

    nrz = int(round((max(cznew)/abs(dz)+0.5)))
    
    Vct3d = []
#    Vct3D = np.zeros((len(cznew), 111, 111))+0.1
    Vct3D = np.zeros(( max(nry), max(nrx), len(cznew)))+0.1
    
#    for i in range(nrz):
    for i in range(np.shape(Vct3D)[0]):

        V0 = V0*dz*(cznew[i]-z0)/dz
 #       VctNew= np.repeat(0.1, nrx[i]*nry[i]).reshape(nry[i], nrx[i])
        R2 = 3.33**2+(6.67-cznew[i])**2
        R = R2**0.5
        j = 0
#        for j in range(nry[i]):
        while j <nry[i]:
            if j <= nryd[i]:
                for k in range(nrxl[i]):
                    Vct3d[j, k] = V0*abs(dx)*(nryd[i]-j)*dy*(nrxl[i]-k)
                    Vct3D[j, k, i] = V0*abs(dx)*(nryd[i]-j)*dy*(nrxl[i]-k)
                for k in range(nrxl[i], nrx[i]):
                    Vct3d[j, k] = V0*abs(dx)*(nryd[i]-j)*dy*(k-nrxl[i])
                    Vct3D[j, k, i] =V0*abs(dx)*(nryd[i]-j)*dy*(k-nrxl[i])
                    if k == 111:
                         break
            else:
                for k in range(nrxl[i]):
                    Vct3D[j, k, i] =V0*abs(dx)*(j-nryd[i])*dy*(nrxl[i]-k)
                    Vct3d[j, k] = V0*abs(dx)*(j-nryd[i])*dy*(nrxl[i]-k)
                for k in range(nrxl[i], nrx[i]):
                    Vct3D[j, k, i] =V0*abs(dx)*(j-nryd[i])*dy*(k-nrxl[i])
                    Vct3d[j, k] = V0*abs(dx)*(j-nryd[i])*dy*(k-nrxl[i])
                    if k == 111:
                         break
            j += 1
#        xnew, ynew = np.mgrid[range(-int(round(nry[i]/2)-0.5), int(round(nry[i]/2)+0.5)), range(-int(round(nry[i]/2)-0.5), int(round(nry[i]/2)+0.5))]
#        x, y = np.mgrid[range(-int(round(nry[i]/2)-0.5), int(round(nry[i]/2)+0.5)), range(-int(round(nrx[i]/2)-0.5), int(round(nrx[i]/2+0.5)+0.5))]
#    tck = scipy.interpolate.griddata(np.meshgrid(range(nry), range(nrx)), Vct, (xnew, ynew), method ='linear')
#        tck = scipy.interpolate.bisplrep(x, y,  Vct3D[i,:,:], s =0)#
#        VctNew = scipy.interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    
        plt.figure()
        plt.pcolor(Vct3D[i, :, :])
        plt.plot(cx[i]+R,cy[i]+R, cx[i]+R,-cy[i]+R)
        plt.plot(3.33+R,6.67+R,'ro')
        plt.title('2D %d th slice extented' % i)
#        Vct3d.append(VctNew)

    xnew = range(max(nrx))
    ynew = range(max(nry))
    znew = range(nrz)
    X,Y = np.meshgrid(xnew, ynew)
    method = 'nearest'
    Vct3Dinterp = scipy.interpolate.griddata(np.mgrid[np.shape(Vct3D)[0], np.shape(Vct3D)[1], np.shape(Vct3D)[2]], Vct3D, (xnew, ynew, znew), method = 'nearest')    
    
    fig = plt.figure()
   # ax = plt.axes(projection ='3d')
    ax = mpl_toolkits.mplot3d.axes3d.Axes3D(fig)
    ax.contour3D(X,Y,VctNew, 50)

#    ax.plot3D(Vct3D[:,0,0:32], Vct3D[:,0:32,0], VctNew[0,0:32,0:32])
    ax.plot3D(Vct3D[0:11,0,0], Vct3D[0,0:11,0], Vct3D[0,0,0:11])
    
    
    uu = []
    xuu = []
    yvv = []
    zuuvv = []
    uu = np.linspace(0, 2*np.pi, int(round(np.pi/(0.14*np.pi)+0.5)))
    for i in range(len(R1)):        
        # x, y, and z are the coordinates of the points for plotting
        # each is arranged in a 8x8 array    
        xuu.append(R1[i]*np.outer(np.cos(uu),np.ones(np.size(uu)))+R1[i])
        yvv.append(R1[i]*np.outer(np.sin(uu),np.ones(np.size(uu)))+R1[i])
        #zuv=10*np.outer(),np.cos(v))
        zuuvv.append(yvv[-1]**2-6.67**2)
    
    xxuu = np.array(xuu)
    yyvv = np.array(yvv) 
    zzuuvv = np.array(zuuvv)  
#    XX, YY, ZZ = np.meshgrid(xxuu[1,:], yyvv[1,:], zzuuvv[1,:])
    
    fig1 = plt.figure()
    #    ax.plot3D(xu,yv,zuv)
    ax1 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig1)
    ax1.plot_wireframe(xxuu[1,:], yyvv[1,:], np.mean(zzuuvv, 2).transpose(),cmap=cm.coolwarm)
#    ax1.contour3D(XX,YY,ZZ)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    
    fig2 = plt.figure()
    #    ax.plot3D(xu,yv,zuv)
    ax2 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig2)
    ax2.plot_surface(xxuu[1,:], yyvv[1,:], np.mean(zzuuvv, 2).transpose(),cmap=cm.coolwarm)
#    ax.plot_surface(np.real(xu),np.real(yv),np.real(zuvext),cmap=cm.coolwarm)
    ax2.plot()
    ax2.set_ylim(-1, 0)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    
    XX, YY =  np.meshgrid(xxuu[1,:], yyvv[1,:])
    XX1 = np.meshgrid(xxuu[1:2,:],xxuu[1:2,:])[0]
    YY1 = np.meshgrid(xxuu[1:2,:],xxuu[1:2,:])[1] 
    ZZ = YY**2-6.67**2
    fig3 = plt.figure()
    ax3 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig3)
#    ax.plot_surface(np.real(xu),np.real(yv),np.real(zuvext),cmap=cm.coolwarm)
    ax3.contour3D(XX,YY,ZZ)
 #   ax.set_ylim(-1, 0)
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('Z')
    
    fig0 = plt.figure()
    ax0 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig0)
#    ax.plot_surface(np.real(xu),np.real(yv),np.real(zuvext),cmap=cm.coolwarm)
    ZZ1 = YY1**2-6.67**2
    ax0.plot3D(np.ravel(YY1),np.ravel(YY),np.ravel(ZZ1))
 #   ax.set_ylim(-1, 0)
    ax0.set_xlabel('X')
    ax0.set_ylabel('Y')
    ax0.set_zlabel('Z')
    
    fig4 = plt.figure()
    ax4 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig4)
    ax4.plot3D(np.ravel(xxuu[1,:]), np.ravel(yyvv[1,:]), np.ravel(np.mean(zzuuvv, 2)))
    ax4.set_ylim(-1, 0)
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_zlabel('Z')
    
    XU = 40*xxuu[1,:]
    YV = 2*abs(40*yyvv[1,:]) -40*yyvv[1,range(7,-1,-1)] 
    YV0 = -20*yyvv[1,range(7,-1,-1)]
    ZUV = 10*np.mean(zzuuvv, 2)+360
    
    fig5 = plt.figure()
    #    ax.plot3D(xu,yv,zuv)
    ax5 = mpl_toolkits.mplot3d.axes3d.Axes3D(fig5)
    ax5.contour3D( 2*abs(ZZ)+20, 40*YY-30, 160*XX,cmap=cm.coolwarm)    
    ax5.contour3D( 2*abs(ZZ)+40, 40*YY-30, 160*XX,cmap=cm.coolwarm)  
    ax5.contour3D( XU, YV, ZUV,cmap=cm.coolwarm)
#    ax.plot_surface(np.real(xu),np.real(yv),np.real(zuvext),cmap=cm.coolwarm)
    ax5.contour3D( XU, YV0, ZUV)
    ax5.plot_surface( yyvv[1,:], 2+8*abs(xxuu[1,:]),np.mean(zzuuvv, 2).transpose())
#    ax5.set_ylim(-1, 0)
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_zlabel('Z')
    
       
    """"
    uu = []
    vv = []
    xuu = []
    yvv = []
    zuuvv = []
    for i in range(round(len(cy[1])/2)): 
        # u is an array from 0 to 2*pi, with 100 elements
        uu.append(np.linspace(0, len(cy[1])*1j, int(round(abs(len(cy[1])*1j)/(2*np.pi)+0.5))))
        # v is an array from 0 to 2*pi, with 100 elements
        vv.append(np.linspace(0, len(cy[1])*1j, int(round(0.5*abs(len(cy[1])*1j)/(np.pi)+0.5))))
        # x, y, and z are the coordinates of the points for plotting
        # each is arranged in a 100x100 array    
        xuu.append((cx[1][0:int(2*np.pi):len(cx[1])]**2+cy[1][0:int(2*np.pi):len(cy[1])]**2)**0.5*np.outer(np.cos(uu),np.sin(vv)))
        yvv.append((cx[1][0:int(2*np.pi):len(cx[1])]**2+cy[1][0:int(2*np.pi):len(cy[1])]**2)**0.5*np.outer(np.sin(uu),np.sin(vv)))
        #zuv=10*np.outer(np.ones((np.size(u))),np.cos(v))
        zuuvv.append(np.outer((cx[1][0:int(2*np.pi):len(cx[1])]**2+cy[1][0:int(2*np.pi):len(cy[1])]**2)**0.5*np.ones(np.size(uu)),np.sin(vv)))
    
     """"

def main():   
        
    f = open("D:\ECG\Person1.txt", "r")
    l = []
    for x in f:
        s = f.readline()
        s =s.split('\t')
        for string in s:
            l.append(string.split())
    f.close()
    
    f1 = open("D:\ECG\Person1_1.txt", "r")
    l1 = []
    for x in f1:
        s = f1.readline()
        s =s.split('\t')
        for string in s:
            l1.append(string.split())
    f1.close()
    
    f2 = open("D:\ECG\Person6.txt", "r")
    l2 = []
    for x in f2:
        s = f2.readline()
        s =s.split('\t')
        for string in s:
            l2.append(string.split())
    f2.close()
    
    f22 = open("D:\ECG\Person6_1.txt", "r")
    l22 = []
    for x in f22:
        s = f22.readline()
        s =s.split('\t')
        for string in s:
            l22.append(string.split())
    f22.close()
    
    lArr=np.array(l[:])
    n_ECG = int(len(lArr)/3)
    ECG = np.repeat(0.1, n_ECG)#filter_data
    ECG_r = np.repeat(0.1, n_ECG)#raw_data
    
    lArr1=np.array(l1[:])
    n_ECG1 = int(len(lArr1)/3)
    ECG1 = np.repeat(0.1, n_ECG1)#filter_data
    ECG_r1 = np.repeat(0.1, n_ECG1)#raw_data
    
    lArr2=np.array(l2[:])
    n_ECG2 = int(len(lArr2)/3)
    ECG2 = np.repeat(0.1, n_ECG2)#filter_data
    ECG_r2 = np.repeat(0.1, n_ECG2)#raw_data
    
    lArr22=np.array(l22[:])
    n_ECG22 = int(len(lArr22)/3)
    ECG22 = np.repeat(0.1, n_ECG22)#filter_data
    ECG_r22 = np.repeat(0.1, n_ECG22)#raw_data
    
    #using zero baseline, filtered ECG
    
    for i in range(2, len(lArr), 3):
        ECG[int(i/3)] = float(lArr[i])  
        ECG_r[int(i/3)] = float(lArr[i-1]) 
    plt.figure()     
    plt.plot(ECG)  
    plt.plot(ECG_r) 
     
        
    for i in range(2, len(l1), 3):
        if str(l1[i])[2] == '-':
            ECG1[int(i/3)] = float(str(l1[i])[2:8])
        else:
            ECG1[int(i/3)] = float(str(l1[i])[2:7])
        if str(l1[i-1])[2] == '-':    
            ECG_r1[int(i/3)] = float(str(l1[i-1])[2:8])  
        else:
            ECG_r1[int(i/3)] = float(str(l1[i-1])[2:7])  
    
    for i in range(2, len(l2), 3):
        if str(l2[i])[2] == '-':
            ECG2[int(i/3)] = float(str(l2[i])[2:8])
        else:
            ECG2[int(i/3)] = float(str(l2[i])[2:7])
        if str(l2[i-1])[2] == '-':    
            ECG_r2[int(i/3)] = float(str(l2[i-1])[2:8])  
        else:
            ECG_r2[int(i/3)] = float(str(l2[i-1])[2:7])  
    
    for i in range(2, len(l22), 3):
        if str(l22[i])[2] == '-':
            ECG22[int(i/3)] = float(str(l22[i])[2:8])
        else:
            ECG22[int(i/3)] = float(str(l22[i])[2:7])
        if str(l22[i-1])[2] == '-':    
            ECG_r22[int(i/3)] = float(str(l22[i-1])[2:8])  
        else:
            ECG_r22[int(i/3)] = float(str(l22[i-1])[2:7])  
    
    
    #    ECG_r1[int(i/3)] = float(lArr1[i-1])
        
    plt.figure()     
    plt.plot(ECG1)  
    plt.plot(ECG_r1) 
    
    fs = 1/20
    fN = int(round((n_ECG/2-1)/2))+1
    fc = 0.5
    #moving windows
    #wk=NT
    RR = 20/13 #s
    nRR = int(np.round(n_ECG/12))
    F = np.fft.fft(ECG)
    plt.plot(F)
    tol = 1
    x = list(ECG[128:int(np.round(128+nRR))])  
    tempy = []
    res = []
    res.append(0)
    tempy1 = []
    res1 = []
    res1.append(0)    
    
    
    for i in range(2*nRR):
        if i == 0:
            temp = np.real(np.fft.ifft(np.fft.fft(0.5*(1-np.cos(2*np.pi/nRR*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))))*np.fft.fft(x)))
            temp1 = np.real(np.fft.ifft(np.fft.fft(0.42-0.5*np.cos(2*np.pi/nRR*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))+0.08*np.cos(4*np.pi/nRR)*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))*np.fft.fft(x))) 
        else:
            T =  int(np.floor(nRR/i/2+1))
    #        if nRR/i/2 <1:
            if nRR/i/4 <1:
                break         
            if len(x) < T:
                while len(x) < T:
                    x.append(0)
        #        temp = np.linspace(-int(np.floor(T/2)), int(np.round(T/2)), T)
            else: 
                for j in range(i): 
                    if j == 0:
                        xtemp = np.zeros((np.size(temp)))
                        xtemp1 = np.zeros((np.size(temp1)))                    
                    xtemp[j*T:(j+1)*T] = np.fft.fft(x[j*T:(j+1)*T])*np.fft.fft(0.5*(1-np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),int(T)))))
                    xtemp1[j*T:(j+1)*T] = np.fft.fft(x[j*T:(j+1)*T])*np.fft.fft(0.42-0.5*np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),int(T))+0.08*np.cos(4*np.pi/T/2/(j+1))*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),T))) 
        #        temp = ''.join([str(np.zeros((len(ECG)-len(temp)))[:]),str(ECG)[:]])     
        #        temp = temp.split() 
    #                temp = np.linspace(min(temp), max(temp),  T) + np.real(np.fft.ifft(np.fft.fft(0.5*(1-np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),T/(j+1)))))*np.fft.fft(xtemp)))        
                    temp += np.real(np.fft.ifft(xtemp)) 
                    temp1 += np.real(np.fft.ifft(xtemp1))                 
                    tempy.append(temp)
                    tempy1.append(temp1)                
                    res.append(abs(np.mean(temp - x)))
                    res1.append(abs(np.mean(temp1 - x)))
                    if abs(np.mean(temp - x)) <= min(res):
                        y = temp
                    if abs(np.mean(temp1 - x)) <= min(res1):
                        y1 = temp1
    plt.figure()
    plt.plot(tempy[-1])
    plt.plot(tempy1[-1])
    
    # Bandpass Filter 
    temp = list(ECG[2:])
    temp.append(ECG[0])
    temp.append(ECG[1])
    temp1 = list(ECG[1:])
    temp1.append(ECG[0])
    
    stage1 = 1/4*(np.array(temp)+2*np.array(temp1)+ECG)
    plt.figure()
    plt.plot(ECG)
    plt.plot(stage1)
    plt.title('stage1 comparison')
    
    Res1 = abs(ECG-stage1)
    SNR1 = []
    SNR1.append(1)
    for i in range(1, np.size(stage1)-1):
        SNR1.append(np.var(stage1[0:i+1])/np.var(Res1[0:i+1]))
    plt.figure()
    plt.plot(Res1)
    plt.title('Res1')
    plt.figure()
    plt.plot(SNR1)
    plt.title('SNR1')
    
    H0 = np.real(np.fft.fft(ECG)[0:int(np.round(0.5*n_ECG))])
    H1 = np.real(np.fft.fft(stage1)[0:int(np.round(0.5*n_ECG))])
    dH = H1/H0
    plt.figure()
    plt.plot(10*np.log(H1)-10*np.log(H0))
    plt.title('Gain')
    
    #Notch Filter
    temp = list(stage1[2:])
    temp.append(stage1[0])
    temp.append(stage1[1])
    temp1 = list(stage1[1:])
    temp1.append(stage1[0])
    
    stage2 = np.array(temp)+2*np.cos(60*np.pi/125)*np.array(temp1)+stage1
    plt.figure()
    plt.plot(ECG)
    plt.plot(stage2)
    plt.title('stage2 comparison')
    
    
    Res2 = abs(ECG-stage2)
    SNR2 = []
    SNR2.append(1)
    for i in range(1, np.size(stage2)-1):
        SNR2.append(np.var(stage2[0:i+1])/np.var(Res2[0:i+1]))
    H2 = np.real(np.fft.fft(stage2)[0:int(np.round(0.5*n_ECG))])
    dH2 = H2/H0
    plt.figure()
    plt.plot(Res2)
    plt.title('Res2')
    plt.figure()
    plt.plot(SNR2)
    plt.title('SNR2')
    plt.figure()
    plt.plot(10*np.log(H2)-10*np.log(H0))
    plt.title('Gain2')
    
    
    #derivative filter
    temp = list(stage2[6:])
    temp.append(stage2[0])
    temp.append(stage2[1])
    temp.append(stage2[2])
    temp.append(stage2[3])
    temp.append(stage2[4])
    temp.append(stage2[5])
    
    stage3 = stage2-np.array(temp)
    plt.figure()
    plt.plot(stage3)
    
    plt.figure()
    plt.plot(ECG)
    plt.plot(stage3)
    plt.title('stage3 comparison')
    
    Res3 = abs(ECG-stage3)
    SNR3 = []
    SNR3.append(1)
    for i in range(np.size(stage3)-1):
        SNR3.append(np.var(stage3[0:i+1])/np.var(Res3[0:i+1]))
    H3 = np.real(np.fft.fft(stage3)[0:int(np.round(0.5*n_ECG))])
    dH3 = H3/H0
    plt.figure()
    plt.plot(Res3)
    plt.title('Res3')
    plt.figure()
    plt.plot(SNR3)
    plt.title('SNR3')
    plt.figure()
    plt.plot(10*np.log(H3)-10*np.log(H0))
    plt.title('Gain3')
    
    
    plt.figure()
    plt.plot(stage3[range(nRR)])
    
    plt.figure()
    plt.plot(stage3[49:65])
    
    QRS = []
    PKID = []
    for i in range(13):
        pks = np.argmax(stage3[49+185*i:195+185*i])
        PKID.append(pks+49+185*i)        
        QRS.append(stage3[PKID[i]-7:PKID[i]+8])
    
    
    Is = 1 
    d = 0.5 
    Ib = np.linspace(1.1, 2.1, 10)#transmembrane current

    t = 10
    V = ECG[128:128+nRR]
    Vt = V[t]
    Pred = Vpredict(Vt, Is, d, Ib, t)
    VPred = Pred[0]
    RohPred = Pred[1]
    iters = 100
    
    source0 = [3.33, 6.67, 6.67, 0, 0, 0]
    dsource = [0.3, 1.9, 2.5, -3.7, 3.2, -34.7]
    #test on t = 10
    reconstructOut = ELBO_simplified(source0, dsource, V, Is, d, t,Ib,iters)
        
    P = reconstructOut[0]
    Q = reconstructOut[1]
    Ib_opt = reconstructOut[2]
    Obj = reconstructOut[3]
    
    plt.figure()
    plt.plot(P)
    plt.plot(Q)
    plt.legend(['based on origin LEAD 1 at N =10','based on approximate J*Ji/Ir/sigma*cos(alpha) at N =10'])
    plt.title('100 simulations at N = 10')

    
    Ppred = []
    Qpred = []
    PVpred = []
    QVpred = []
    Opt = []
    IbOpt = []

    #reconstruct magnitude
    for t in range(nRR):
        reconstructOut = ELBO_simplified(source0, dsource, V, Is, d, t,Ib,iters)            
        P = reconstructOut[0]
        Q = reconstructOut[1]
        Ib_opt = reconstructOut[2]
        Obj = reconstructOut[3]
        Ppred.append(np.mean(P))
        Qpred.append(np.mean(Q))
        if Ppred[-1]*V[t]<0:
            PVpred.append(-Ppred[-1])    
        else:
            PVpred.append(Ppred[-1]) 
        if Qpred[-1]*V[t]<0:
            QVpred.append(-Qpred[-1])
        else:
            QVpred.append(Qpred[-1])     
        IbOpt.append(Ib_opt)
        Opt.append(Obj)    

    #sensitivity analysis        
    plt.figure()
    plt.plot(abs(V))
    plt.plot(Ppred)
    plt.plot(Qpred)
    plt.legend(['Origin','based on origin LEAD 1 at first QRS','based on approximate J*Ji/Ir/sigma*cos(alpha) at first QRS'])
    plt.title('100 simulations mean magnitude at first QRS')

    plt.figure()
    plt.plot(V)
    plt.plot(PVpred)
    plt.plot(QVpred)
    plt.legend(['OriginVoltage','Voltage based on origin LEAD 1 at first QRS','Voltage based on approximate J*Ji/Ir/sigma*cos(alpha) at first QRS'])
    plt.title('100 simulations mean Voltage at first QRS')
    
    ResP = PVpred - V
    ResQ = QVpred - V
    
    Pmean = []
    Qmean = []
    PCV = []
    QCV = []
    PResMean = []
    QResMean = []
    PResCV = []
    QResCV = [] 

    for t in range(nRR-1):    
        Pmean.append(np.mean(PVpred[0:(t+1)]))
        Qmean.append(np.mean(QVpred[0:(t+1)]))
        PResMean.append(np.mean(ResP[0:(t+1)]))
        QResMean.append(np.mean(ResQ[0:(t+1)]))        
        if t == 0:
            PCV.append(0)
            QCV.append(0)
            PResCV.append(0)
            QResCV.append(0)
        else:
            PCV.append(np.std(PVpred[0:(t+1)])/(1+Pmean[-1]))
            QCV.append(np.std(QVpred[0:(t+1)])/(1+Qmean[-1]))
            PResCV.append(np.std(ResP[0:(t+1)])/(1+PResMean[-1]))
            QResCV.append(np.std(ResQ[0:(t+1)])/(1+QResMean[-1])) 
    
    Pmean.append(Pmean[-1])
    Qmean.append(Qmean[-1])
    PCV.append(PCV[-1])
    QCV.append(QCV[-1])
    PResMean.append(PResMean[-1])
    QResMean.append(QResMean[-1])
    PResCV.append(PResCV[-1])
    QResCV.append(QResCV[-1])
    
    
    plt.figure()
    plt.plot(ResP) 
    x = range(nRR)
    y = PResMean
    yerr = PResCV
    plt.errorbar(x, y, yerr = yerr[::-1])
    plt.legend(['Res based on origin LEAD 1 at first QRS'])
    plt.title('100 simulations Residule mean and variational at first QRS')
    
    plt.figure()
    plt.plot(ResQ)
    y = QResMean
    yerr = QResCV
    plt.errorbar(x, y, yerr = yerr[::-1])    
    plt.legend(['Res based on approximate J*Ji/Ir/sigma*cos(alpha) at first QRS'])
    plt.title('100 simulations Residule mean and variational at first QRS')

    plt.figure()
    y = Pmean
    yerr = PCV
    plt.plot(Pmean)   
    plt.errorbar(x, y, yerr = yerr[::-1])
    plt.legend(['Res based on origin LEAD 1 at first QRS'])
    plt.title('100 simulations Residule mean and variational at first QRS')

    plt.figure()
    y = Qmean
    yerr = QCV
    plt.plot(Qmean)
    plt.errorbar(x, y, yerr = yerr[::-1])    
    plt.legend(['Mean based on approximate J*Ji/Ir/sigma*cos(alpha) at first QRS'])
    plt.title('100 simulations Predict mean and variational at first QRS')

    

    
    plt.figure()
    plt.plot(np.transpose(np.array(QRS)))
    plt.legend()
    plt.title('QRS')
    
    QRSdf = pd.DataFrame(QRS)
 #   rp.summary_cont(QRSdf['stats'])
    RP = rp.summary_cont(QRSdf.transpose())
    QRS_ANOVA = sci.f_oneway(RP.iloc[0:4,2:7], RP.iloc[4:8,2:7],  RP.iloc[8:13,2:7])
    
    lm = LinearRegression()
    lm = sfa.ols('outcome_variable ~ C(group_variable)', data = [RP.iloc[0:4,2:7], RP.iloc[4:8,2:7], RP.iloc[8:13,2:7]]).fit()   
    
    
    #Main detector
    #The event analysis is done if there is at least one
    #crossing followed by 180 ms(2500/20*180/60=375 samples) without a new crossing or if
    #there are more than 4 crossings.
    crs = 0
    BT = np.mean(stage3) + 0.12
    DT = np.mean(stage3) - 0.01
    nQRS = 0
    nNoise = 0
    tempidx = []
    QRS = []
    yE = []
    PTs = []
    
    for idx in range(len(stage3)):
    #    if stage3[idx] < DT: flag = -1
    #    elif stage3[idx] > BT: flag = 1
    #    else: flag = 0 
    
        if (stage3[idx] - DT)*(DT - stage3[idx+1]) >= 0 + (stage3[idx] - BT)*(BT - stage3[idx+1]) >= 0:
            tempidx.append((stage3[idx] - DT)*(DT - stage3[idx+1]) >= 0 + (stage3[idx] - BT)*(BT - stage3[idx+1]) >= 0)
            crs += 1
            PK = max(stage3[(idx+1):(idx+375)])
            PKID = np.argmax(stage3[(idx+1):(idx+375)])
            Pts = CrossingCount(stage3[(idx+1):(idx+375)], idx, BT, DT, nQRS, nNoise, PK, PKID)
            if Pts.name == 'QRS':
                QRS.append(Pts)
                nQRS += 1
                BT = Pts.BT
                DT = Pts.DT
            elif Pts.name == 'noise':
                nNoise += 1
            elif Pts.name == 'DK':
                if np.size(yE)==0:
                    E = stage3**2
                    for i in range(1, 20):
                        temp = list(stage3[i:])
                        while j in range(i):
                            temp.append(j)
                        E += temp**2
                    yE = E
                if EnergyCrossingCount(E[(idx+1):(idx+375)], Pts, nQRS, nNoise, PK, PKID) == False:
                    Pts.name = 'Else'
            PTs.append(Pts) 
            
    
    
if __name__ == "__main__":
    main()       
                


##

for i in range(2*nRR):
    if i == 0:
        temp = np.real(np.fft.ifft(np.fft.fft(0.5*(1-np.cos(2*np.pi/nRR*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))))*np.fft.fft(x)))
        temp1 = np.real(np.fft.ifft(np.fft.fft(0.42-0.5*np.cos(2*np.pi/nRR*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))+0.08*np.cos(4*np.pi/nRR)*np.linspace(-int(np.floor(nRR/2)),int(np.floor(nRR/2)),nRR))*np.fft.fft(x))) 
    else:
        T =  int(np.floor(nRR/i/2+1))
#        if nRR/i/2 <1:
        if nRR/i/4 <1:
            break         
        if len(x) < T:
            while len(x) < T:
                x.append(0)
    #        temp = np.linspace(-int(np.floor(T/2)), int(np.round(T/2)), T)
        else: 
            for j in range(i): 
                if j == 0:
                    xtemp = np.zeros((np.size(temp)))
                    xtemp1 = np.zeros((np.size(temp1)))                    
                xtemp[j*T:(j+1)*T] = np.fft.fft(x[j*T:(j+1)*T])*np.fft.fft(0.5*(1-np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),int(T)))))
                xtemp1[j*T:(j+1)*T] = np.fft.fft(x[j*T:(j+1)*T])*np.fft.fft(0.42-0.5*np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),int(T))+0.08*np.cos(4*np.pi/T/2/(j+1))*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),T))) 
    #        temp = ''.join([str(np.zeros((len(ECG)-len(temp)))[:]),str(ECG)[:]])     
    #        temp = temp.split() 
#                temp = np.linspace(min(temp), max(temp),  T) + np.real(np.fft.ifft(np.fft.fft(0.5*(1-np.cos(2*np.pi/T*np.linspace(-int(np.floor(T/2/(j+1))),int(np.floor(T/2/(j+1))),T/(j+1)))))*np.fft.fft(xtemp)))        
                temp += np.real(np.fft.ifft(xtemp)) 
                temp1 += np.real(np.fft.ifft(xtemp1))                 
                tempy.append(temp)
                tempy1.append(temp1)                
                res.append(abs(np.mean(temp - x)))
                res1.append(abs(np.mean(temp1 - x)))
                if abs(np.mean(temp - x)) <= min(res):
                    y = temp
                if abs(np.mean(temp1 - x)) <= min(res1):
                    y1 = temp1


#        0.5*(1-np.cos(2*np.pi/M*np.linspace(-M/2-1,M/2+1,M+4)*np.exp(-2j*pi*j**2/M)))
#plt.plot(np.fft.fft(0.5*(1-np.cos(2*np.pi/M*np.linspace(-M/2,M/2,M)))))
#plt.plot(np.fft.ifft(np.fft.fft(0.5*(1-np.cos(2*np.pi*(M/j/2-M/j)/(M/j)*np.linspace(-int(np.floor(M/j+2)),int(np.floor(M/j+2)),int(np.floor(M/j+2))))))))
#plt.plot(np.fft.fft(0.5*(1-np.cos(2*np.pi*(M/j/2-M/j)/(M/j)*np.linspace(-int(np.floor(M/j+2)),int(np.floor(M/j+2)),int(np.floor(M/j+2)))))))
#plt.plot(np.fft.fft(ECG[range(int(np.floor(M/j+2)))]))  
#plt.plot(0.5*(1-np.cos(2*np.pi/M*np.linspace(-M/2,M/2,M))))
plt.plot(x)
plt.figure()
plt.plot(res)
plt.figure()
plt.plot(y)
plt.figure()
plt.plot(res1)
plt.figure()
plt.plot(y1)
    
 
#    temp += ECG[] 
#    H.append(temp) 
#    if np.real(np.mean(ECG[range(2*i+1)] - temp))< tol:
#        break
#    if np.mod(i,500)==0:
#        plt.figure()
#       plt.plot(temp)

 
#Add drift baseline b(t) = C*sum(a_k*cos(2*pi*df*t)+phi_k)
