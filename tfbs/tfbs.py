import pandas as pd
pd.options.display.float_format = '${:,.8f}'.format
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import unicodedata
import os
import scipy.special as sps

def prune(P0, j, L, Seq, k, LDNA, l, c, fb):
    Pf_ = np.array(list([0]*len(x)))
    Pb_ = np.array(list([0]*len(x)))
    L = L-1;
    if fb == 'f' :
        if k - 1 ==  k - L:
            nfn_ = Seq[k-1].count(LDNA[j])
        else:
            nfn_ = Seq[k-1:k-L+1].count(LDNA[j])
        if nfn_ == 0:
            nfn_ = Seq[k:k+2*L-1+1+1].count(LDNA[j])/len(LDNA[j])**(l+2)
        if k ==  k - L:
            nfd_ = Seq[k].count(LDNA[j])
        else:
            nfd_ = Seq[k-L:k+1].count(LDNA[j])
        if nfd_ == 0:
            nfd_ = Seq[k:k+2*L-1+1].count(LDNA[j])/len(LDNA[j])**(l+2)
        if nfd_ == 0:
            if nfn_ == 0:
                Pf_ = 0
            else:
                Pf_ = 0.25
        else:
            Pf_ = nfn_/nfd_   
        
        KL = sum(Pf_*np.log(2, P0/Pf))   
        if KL <= c:
            Pf_= 0
        P = Pf_
            
        
    else:
        if k - 1 ==  k - Lb:
            nbn_ = Seq[k-1].count(dDNA[j])
        else:
            nbn_ = Seq[k-L:k-1+1].count(dDNA[j])
        if nbn_ == 0:
            nbn_ = Seq[k:k+2*L-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)
        if k ==  k - L:
            nbd_ = Seq[k-1].count(dDNA[j])
        else:
            nbd_ = Seq[k-L:k+1].count(dDNA[j])
        if nbd_ == 0:
            nbd_ = Seq[k:k+2*L-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)
        if nbd_ == 0:
            if nbn_ == 0:
                Pb_ = 0
            else:
                Pb_ = 0.25
        else:
            Pb_ = nbn_/nbd_   
        
        KL = sum(Pb_*np.log(2, P0/Pb))  
        psi = len(LDNA)**(l+1)/nbn_
        while (L > 1) & (KL <= c*psi):
            L = L - 1
            Pb_ = prune(Pb_, j, L, Seq, k, dDNA, l, c, 'b')
        P = Pb_
    return P
#path = r'D:/Introduction to R/HMM//'
file = r'P:/GCMC/GCMC/STXBP1pheno.txt'
#file = r'D:/Introduction to R/HMM/STXBP1pheno.txt'
dataset = pd.read_csv(file, sep = ',')
Seq = str(dataset.x[1])
seq = list(dataset.x[1])
x = seq[0:99] #example 100
        

dDNA = ['A','T','C','G']
Lf = 1
Lb = 5# binding site
c = 0.65


Pf = np.array(list([0.0]*len(dDNA)))
Pb = np.array(list([0.0]*len(dDNA)))
nfn = np.array(list([0.0]*len(dDNA)))
nfd = np.array(list([0.0]*len(dDNA)))
nbn = np.array(list([0.0]*len(dDNA)))
nbd = np.array(list([0.0]*len(dDNA)))
PF = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
PB = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))

pred_f =np.array(list(['']*len(x)))
accAf = np.array(list([0.0]*len(x)))
accCf = np.array(list([0.0]*len(x)))
accTf = np.array(list([0.0]*len(x)))
accGf = np.array(list([0.0]*len(x)))
accf = np.array([accAf,accCf,accTf,accGf])
tfbsAf = np.array(list([0.0]*len(x)))
tfbsTf = np.array(list([0.0]*len(x)))
tfbsCf = np.array(list([0.0]*len(x)))
tfbsGf = np.array(list([0.0]*len(x))) 


#L = 1 find tfbs with 10 box(fragment with significantly higher accuracy)
for i in np.arange(1, len(x)):  
    for j in range(len(dDNA)):
        if i - 1 ==  i - Lf:
            nfn[j] = Seq[i-1].count(dDNA[j])
        else:
            #[]: upper bound excluded
            nfn[j] = Seq[i-Lf:i-1+1].count(dDNA[j])
        if nfn[j] == 0:
            nfn[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1) 
        if i ==  i - Lf:
            nfd[j] = Seq[i].count(dDNA[j]) 
        else:
            nfd[j] = Seq[i-Lf:i+1].count(dDNA[j])
        if nfd[j] == 0:
            nfd[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)
        if nfd[j] == 0:
            if nfn[j] == 0:
                Pf[j] = 0
            else:
                Pf[j] = 0.25
        else:
            Pf[j] = nfn[j]/nfd[j]        
        if  x[i] == dDNA[0]:
            PF[0,j] = Pf[j]
    #            PB[0,j] = max(Pb[j], Pb_[j])               
        if  x[i] == dDNA[1]: 
            PF[1,j] = Pf[j]
    #            PB[1,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[2]: 
            PF[2,j] = Pf[j]
    #            PB[2,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[3]: 
            PF[3,j] = Pf[j]
    #            PB[3,j] = max(Pb[j], Pb_[j])     
        if j == 3:
            nfn = np.array(list([0]*len(dDNA)))
            nfd = np.array(list([0]*len(dDNA)))
    P0 = np.dot([1,1,1,1], PF)
    pred_f[i] = dDNA[np.argmax(P0)]
    accAf[i] = sum((pred_f[0:i+1]==x[0:i+1]) * (pred_f[0:i+1] == np.repeat(dDNA[0],i+1)))/(i+1)* len(dDNA)
    accTf[i] = sum((pred_f[0:i+1]==x[0:i+1]) * (pred_f[0:i+1] == np.repeat(dDNA[1],i+1)))/(i+1)* len(dDNA)
    accCf[i] = sum((pred_f[0:i+1]==x[0:i+1]) * (pred_f[0:i+1] == np.repeat(dDNA[2],i+1)))/(i+1)* len(dDNA)
    accGf[i] = sum((pred_f[0:i+1]==x[0:i+1]) * (pred_f[0:i+1] == np.repeat(dDNA[3],i+1)))/(i+1)* len(dDNA)        
accAf[accAf>1] = 1
accTf[accTf>1] = 1
accCf[accCf>1] = 1
accGf[accGf>1] = 1
accf = np.array([accAf,accCf,accTf,accGf])
tfbsAf[(accAf>np.mean(accAf))*(accAf>=accCf)*(accAf>=accTf)*(accAf>=accGf)*(pred_f==dDNA[0])] = 1
tfbsTf[(accTf>np.mean(accTf))*(accTf>=accAf)*(accTf>=accCf)*(accTf>=accGf)*(pred_f==dDNA[1])] = 1
tfbsCf[(accCf>np.mean(accCf))*(accCf>=accAf)*(accCf>=accTf)*(accCf>=accGf)*(pred_f==dDNA[2])] = 1
tfbsGf[(accGf>np.mean(accGf))*(accGf>=accAf)*(accGf>=accCf)*(accGf>=accTf)*(pred_f==dDNA[3])] = 1
tfbsf = np.array([tfbsAf,tfbsCf,tfbsTf,tfbsGf]) 
##L = 5 find tfbs with 10 box(fragment with significantly higher accuracy)
Pf = np.array(list([0.0]*len(dDNA)))
Pb = np.array(list([0.0]*len(dDNA)))
nfn = np.array(list([0.0]*len(dDNA)))
nfd = np.array(list([0.0]*len(dDNA)))
nbn = np.array(list([0.0]*len(dDNA)))
nbd = np.array(list([0.0]*len(dDNA)))
PF = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
PB = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))

pred_b =np.array(list(['']*len(x)))
accAb = np.array(list([0.0]*len(x)))
accCb = np.array(list([0.0]*len(x)))
accTb = np.array(list([0.0]*len(x)))
accGb = np.array(list([0.0]*len(x)))
accb = np.array([accAb,accCb,accTb,accGb])
tfbsAb = np.array(list([0.0]*len(x)))
tfbsTb = np.array(list([0.0]*len(x)))
tfbsCb = np.array(list([0.0]*len(x)))
tfbsGb = np.array(list([0.0]*len(x))) 

for i in np.arange(Lb, len(x)):  
    for j in range(len(dDNA)):        
        if i - 1 ==  i - Lb:
            nbn[j] = Seq[i-1].count(dDNA[j])
        else:
            #[]: upper bound excluded
            nbn[j] = Seq[i-Lb:i-1+1].count(dDNA[j])
        if nbn[j] == 0:
            nbn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(i+1) 
        if i ==  i - Lb:
            nbd[j] = Seq[i].count(dDNA[j]) 
        else:
            nbd[j] = Seq[i-Lb:i+1].count(dDNA[j])
        if nbd[j] == 0:
            nbd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)
        if nbd[j] == 0:
            if nbn[j] == 0:
                Pb[j] = 0
            else:
                Pb[j] = 0.25
        else:
            Pb[j] = nbn[j]/nbd[j]        
        if  x[i] == dDNA[0]:
            PB[0,j] = Pb[j]
    #            PB[0,j] = max(Pb[j], Pb_[j])               
        if  x[i] == dDNA[1]: 
            PB[1,j] = Pb[j]
    #            PB[1,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[2]: 
            PB[2,j] = Pb[j]
    #            PB[2,j] = max(Pb[j], Pb_[j])
        if  x[i] == dDNA[3]: 
            PB[3,j] = Pb[j]
    #            PB[3,j] = max(Pb[j], Pb_[j])     
        if j == 3:
            nbn = np.array(list([0]*len(dDNA)))
            nbd = np.array(list([0]*len(dDNA)))    
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
    P0 = np.dot([1,1,1,1], PB)
    pred_b[i] = dDNA[np.argmax(P0)]
    accAb[i] = sum((pred_b[0:i+1]==x[0:i+1]) * (pred_b[0:i+1] == np.repeat(dDNA[0],i+1)))/(i+1)* len(dDNA)
    accTb[i] = sum((pred_b[0:i+1]==x[0:i+1]) * (pred_b[0:i+1] == np.repeat(dDNA[1],i+1)))/(i+1)* len(dDNA)
    accCb[i] = sum((pred_b[0:i+1]==x[0:i+1]) * (pred_b[0:i+1] == np.repeat(dDNA[2],i+1)))/(i+1)* len(dDNA)
    accGb[i] = sum((pred_b[0:i+1]==x[0:i+1]) * (pred_b[0:i+1] == np.repeat(dDNA[3],i+1)))/(i+1)* len(dDNA)     

accb = np.array([accAb,accCb,accTb,accGb])
tfbsAb[(accAb>np.mean(accAb))*(accAb>=accCb)*(accAb>=accTb)*(accAb>=accGb)*(pred_b==dDNA[0])] = 1
tfbsTb[(accTb>np.mean(accTb))*(accTb>=accAb)*(accTb>=accCb)*(accTb>=accGb)*(pred_b==dDNA[1])] = 1
tfbsCb[(accCb>np.mean(accCb))*(accCb>=accAb)*(accCb>=accTb)*(accCb>=accGb)*(pred_b==dDNA[2])] = 1
tfbsGb[(accGb>np.mean(accGb))*(accGb>=accAb)*(accGb>=accCb)*(accGb>=accTb)*(pred_b==dDNA[3])] = 1
tfbsb = np.array([tfbsAb,tfbsCb,tfbsTb,tfbsGb]) 


plt.figure
plt.subplot(221)
plt.plot(accf.T)
plt.title('acc with L = 1')
plt.legend(['A','T','C','G'])
plt.subplot(222)
plt.plot(tfbsf.T)
plt.title('accordant tfbs with L = 1')
plt.legend(['A','T','C','G'])
plt.subplot(223)
plt.plot(accb.T)
plt.title('acc with L = 5')
plt.legend(['A','T','C','G'])
plt.subplot(224)
plt.plot(tfbsb.T)
plt.title('accordant tdbs with L = 1')
plt.legend(['A','T','C','G'])

Pf = np.array(list([0.0]*len(dDNA)))
Pb = np.array(list([0.0]*len(dDNA)))
nfn = np.array(list([0.0]*len(dDNA)))
nfd = np.array(list([0.0]*len(dDNA)))
nbn = np.array(list([0.0]*len(dDNA)))
nbd = np.array(list([0.0]*len(dDNA)))
PF = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
PB = np.repeat(np.zeros(len(dDNA)), len(dDNA)).reshape(len(dDNA), len(dDNA))
pred_ =np.array(['']*(len(x)))
t = list([[]]*(len(x)))
parents = list([[]]*(len(x)))
childsN =list([[]]*(len(x)))

tfbs = np.array(list([0.0]*np.size(tfbsb,1)))
tfbs[sum(tfbsf == np.zeros([4,99]),0)!=4] = 1

#fix L, 10 datapoints, homogenious VOM(5,0.65) tree
for i in np.arange(1, len(x)-1):
    if tfbs[i] == 0 | i < Lb:
        for j in range(len(dDNA)):
            if i - 1 ==  i - Lf:
                nfn[j] = Seq[i-1].count(dDNA[j])
            else:
                nfn[j] = Seq[i-Lf:i-1+1].count(dDNA[j])
            if nfn[j] == 0:
                nfn[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)
            if i ==  i - Lf:
                nfd[j] = Seq[i].count(dDNA[j])
            else:       
                nfd[j] = Seq[i-Lf:i+1].count(dDNA[j])
            if nfd[j] == 0:
                nfd[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1)
            if nfd[j] == 0:
                if nfn[j] == 0:
                    Pf[j] = 0
                else:
                    Pf[j] = 0.25
            else:
                Pf[j] = nfn[j]/nfd[j]            
            if  x[i] == dDNA[0]:
                PF[0,j] = Pf[j]
            #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PF[1,j] = Pf[j]
            #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PF[2,j] = Pf[j]
            #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PF[3,j] = Pf[j] 
            #            PB[3,j] = max(Pb[j], Pb_[j])  
            if j == 3:
                nfn = np.array(list([0]*len(dDNA)))
                nfd = np.array(list([0]*len(dDNA)))
    #            nbn_ = np.array(list([0]*len(LDNA)))
    #            nbd_ = np.array(list([0]*len(LDNA)))
    #            nbn = np.array(list([0]*len(LDNA)))
    #            nbd = np.array(list([0]*len(LDNA)))             
    #    t.append(np.dot(P, [1,1,1,1])) 
        P0 = np.dot([1,1,1,1], PF/len(dDNA))
        t[i].append(P0)
        childsN[i].append(4)
        for j in range(len(dDNA)):    
            if i - 1 ==  i - Lf +1:
                nfn[j] = Seq[i-1].count(dDNA[j])
            else:
                nfn[j] = Seq[i-Lf+1:i-1+1].count(dDNA[j])
            if nfn[j] == 0:
                nfn[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1) 
            if i ==  i - Lf:
                nfd[j] = Seq[i].count(dDNA[j])
            else:       
                nfd[j] = Seq[i-Lf+1:i+1].count(dDNA[j])
            if nfd[j] == 0:
                nfd[j] = Seq[i:i+2*Lf-1+1].count(dDNA[j])/len(dDNA[j])**(i+1) 
            if nfd[j] == 0:
                if nfn[j] == 0:
                    Pf[j] = 0
                else:
                    Pf[j] = 0.25
            else:
                Pf[j] = nfn[j]/nfd[j]           
            KL = sum(Pf*np.log(2, P0/Pf))
            if KL <= c:
                Pf[j] = prune(Pf[j], j, Lf, Seq, i, dDNA, 1, c, 'f')        
            if x[i] == dDNA[0]:
                PF[0,j] = Pf[j]
            #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PF[1,j] = Pf[j]
            #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PF[2,j] = Pf[j]
            #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PF[3,j] = Pf[j]
            #            PB[3,j] = max(Pb[j], Pb_[j])     
            if j == 3:
                nfn = np.array(list([0]*len(dDNA)))
                nfd = np.array(list([0]*len(dDNA)))
        temp = np.dot([1,1,1,1], PB)
        ind = np.array([0,1,2,3])
#                parents.append(dDNA[ind*(temp!= np.array([0.0]*len(dDNA)))])
                
        for j in ind*(temp!= np.array([0.0]*len(dDNA))): 
            if sum(PF[0:,j] >np.ones(4))>0:
                t[i].append(PF[0:,j]/(max(PF[0:,j])*len(dDNA)))
            parents[i].append(dDNA[j])
            childsN[i].append(sum(temp!= 0.0)) 
        pred_[i] = dDNA[np.argmax(sum(np.array(t[i]),1))]

            
    else:
        l = 0
        t[i] = list([[]]*len(dDNA))
        parents[i] = list([[]]*len(dDNA))
        childsN[i] =list([[]]*len(dDNA))
        for j in range(len(dDNA)):
            if i - 1 ==  i - Lb:
                nbn[j] = Seq[i-1].count(dDNA[j])
            else:
                nbn[j] = Seq[i-Lb:i-1+1].count(dDNA[j])
            if nbn[j] == 0:
                nbn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)
            if i  ==  i - Lf:
                nbd[j] = Seq[i].count(dDNA[j])
            else:
                nbd[j] = Seq[i-Lb:i+1].count(dDNA[j])
            if nbd[j] == 0:
                nbd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2) 
            if nbd[j] == 0:
                if nbn[j] == 0:
                    Pf[j] = 0
                else:
                    Pf[j] = 0.25
            else:
                Pb[j] = nbn[j]/nbd[j]            
            if  x[i] == dDNA[0]:
                PB[0,j] = Pb[j]
    #            PB[0,j] = max(Pb[j], Pb_[j])               
            if  x[i] == dDNA[1]: 
                PB[1,j] = Pb[j]
    #            PB[1,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[2]: 
                PB[2,j] = Pb[j]
    #            PB[2,j] = max(Pb[j], Pb_[j])
            if  x[i] == dDNA[3]: 
                PB[3,j] = Pb[j]
    #            PB[3,j] = max(Pb[j], Pb_[j])     
            if j == 3:
                nbn = np.array(list([0]*len(dDNA)))
                nbd = np.array(list([0]*len(dDNA)))
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
        P0 = np.dot([1,1,1,1], PB)
        if l == 0:
            t[i].append(P0/len(dDNA))
            childsN[i].append(4)
            l += 1
            for j in range(len(dDNA)):
                parents[i].append(dDNA[j])
                if sum(PB[0:,j] >np.ones(4))>0:
                    t[i].append(PB[0:,j]/(max(PB[0:,j])*len(dDNA)))
                childsN[i].append(4)        
                for j in range(len(dDNA)):
                    if i - 1 ==  i - Lb +1 :
                        nbn[j] = Seq[i-1].count(dDNA[j])
                    else:
                        nbn[j] = Seq[i-Lb+1:i-1+1].count(dDNA[j])
                    if nbn[j] == 0:
                        nbn[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)
                    if i  ==  i - Lf:
                        nbd[j] = Seq[i].count(dDNA[j])
                    else:
                        nbd[j] = Seq[i-Lb+1:i+1].count(dDNA[j])
                    if nbd[j] == 0:
                        nbd[j] = Seq[i:i+2*Lb-1+1].count(dDNA[j])/len(dDNA[j])**(l+2)
                    if nbd[j] == 0:
                        if nbn[j] == 0:
                            Pb[j] = 0
                        else:
                            Pb[j] = 0.25
                    else:
                        Pb[j] = nbn[j]/nbd[j]  
                
                    KL = sum(Pb*np.log(2, P0/Pb)) 
                    psi = len(dDNA)**(l+1)/nbn[j]
                    if KL <= c*psi:
                        Pb[j] = prune(Pb[j], j, Lb, Seq, i, dDNA, l, c, 'b')           
                    if x[i] == dDNA[0]:
                        PB[0,j] = Pb[j]
                #            PB[0,j] = max(Pb[j], Pb_[j])               
                    if  x[i] == dDNA[1]: 
                        PB[1,j] = Pb[j]
                #            PB[1,j] = max(Pb[j], Pb_[j])
                    if  x[i] == dDNA[2]: 
                        PB[2,j] = Pb[j]
                #            PB[2,j] = max(Pb[j], Pb_[j])
                    if  x[i] == dDNA[3]: 
                        PB[3,j] = Pb[j]
                #            PB[3,j] = max(Pb[j], Pb_[j])     
                    if j == 3:
                        nbn = np.array(list([0]*len(dDNA)))
                        nbd = np.array(list([0]*len(dDNA)))
                temp = np.dot([1,1,1,1], PB)
                l += 1
                ind = np.array([0,1,2,3])
#                parents.append(dDNA[ind*(temp!= np.array([0.0]*len(dDNA)))])
                
                for j in ind*(temp!= np.array([0.0]*len(dDNA))):  
                    if sum(PB[0:,j] >np.ones(4))>0:
                        t[i].append(PB[0:,j]/(max(PB[0:,j])*len(dDNA)))
                    parents[i].append(dDNA[j])
                childsN[i].append(sum(temp!= 0.0))
                temp = t[i]
                for tp in range(len(temp)-1):
                    if temp[tp] == []:
                        temp[tp] = np.zeros(len(dDNA))
                t[i] = temp
                pred_[i] = dDNA[np.argmax(sum(t[i],1))]
                
pred_[len(x)-1] = x[len(x)-1]
ACCb = sum(pred_b ==x[0:(len(pred_b))])/len(pred_b)
meanACCb = np.zeros(len(pred_b))
stdACCb = np.zeros(len(pred_b))
ACCf = sum(pred_f ==x[0:(len(pred_f))])/len(pred_f) 
meanACCf = np.zeros(len(pred_f))
stdACCf = np.zeros(len(pred_f))    
ACC = sum(pred_ ==x[0:(len(pred_))])/len(pred_)
meanACC = np.zeros(len(pred_))
stdACC = np.zeros(len(pred_))
        
accA_ = np.array(list([0.0]*len(pred_)))
accC_ = np.array(list([0.0]*len(pred_)))
accT_ = np.array(list([0.0]*len(pred_)))
accG_ = np.array(list([0.0]*len(pred_)))
TPA_ = np.array(list([0.0]*len(pred_)))
TPC_ = np.array(list([0.0]*len(pred_)))
TPT_ = np.array(list([0.0]*len(pred_)))
TPG_ = np.array(list([0.0]*len(pred_)))
TNA_ = np.array(list([0.0]*len(pred_)))
TNC_ = np.array(list([0.0]*len(pred_)))
TNT_ = np.array(list([0.0]*len(pred_)))
TNG_ = np.array(list([0.0]*len(pred_)))
acc_ = np.array([accAf,accCf,accTf,accGf])
tfbsA_ = np.array(list([0.0]*len(pred_)))
tfbsT_ = np.array(list([0.0]*len(pred_)))
tfbsC_ = np.array(list([0.0]*len(pred_)))
tfbsG_ = np.array(list([0.0]*len(pred_))) 
for i in np.arange(0, len(pred_)-1):  
    accA_[i] = sum((pred_[0:i+1]==x[0:i+1]) * (pred_[0:i+1] == np.repeat(dDNA[0],i+1)))/(i+1)* len(dDNA)
    accT_[i] = sum((pred_[0:i+1]==x[0:i+1]) * (pred_[0:i+1] == np.repeat(dDNA[1],i+1)))/(i+1)* len(dDNA)
    accC_[i] = sum((pred_[0:i+1]==x[0:i+1]) * (pred_[0:i+1] == np.repeat(dDNA[2],i+1)))/(i+1)* len(dDNA)
    accG_[i] = sum((pred_[0:i+1]==x[0:i+1]) * (pred_[0:i+1] == np.repeat(dDNA[3],i+1)))/(i+1)* len(dDNA)     
    TPA_[i] = sum( (pred_[0:i+1] == np.repeat(dDNA[0],i+1)))/sum(x[0:i+1] == np.repeat(dDNA[0],i+1))
    TPT_[i] = sum((pred_[0:i+1] == np.repeat(dDNA[1],i+1)))/sum(x[0:i+1] == np.repeat(dDNA[1],i+1))
    TPC_[i] = sum((pred_[0:i+1] == np.repeat(dDNA[2],i+1)))/sum(x[0:i+1] == np.repeat(dDNA[2],i+1))
    TPG_[i] = sum((pred_[0:i+1] == np.repeat(dDNA[3],i+1)))/sum(x[0:i+1] == np.repeat(dDNA[3],i+1))    
    TNA_[i] = sum((pred_[0:i+1] != np.repeat(dDNA[0],i+1)))/sum(x[0:i+1] != np.repeat(dDNA[0],i+1))
    TNT_[i] = sum((pred_[0:i+1] != np.repeat(dDNA[1],i+1)))/sum(x[0:i+1] != np.repeat(dDNA[1],i+1))
    TNC_[i] = sum((pred_[0:i+1] != np.repeat(dDNA[2],i+1)))/sum(x[0:i+1] != np.repeat(dDNA[2],i+1))
    TNG_[i] = sum((pred_[0:i+1] != np.repeat(dDNA[3],i+1)))/sum(x[0:i+1] != np.repeat(dDNA[3],i+1))     


accA_[accA_>1] = 1
accT_[accT_>1] = 1
accC_[accC_>1] = 1
accG_[accG_>1] = 1

acc_ = np.array([accA_,accC_,accT_,accG_])
TP_ = np.array([TPA_/max(TPA_),TPC_/max(TPC_),TPT_/max(TPT_),TPG_/max(TPG_)])
TN_ = np.array([TNA_/max(TNA_),TNC_/max(TPC_),TNT_/max(TNT_),TNG_/max(TNG_)])

tfbsA_[(accA_>np.mean(accA_))*(accA_>=accC_)*(accA_>=accT_)*(accA_>=accG_)*(pred_==dDNA[0])] = 1
tfbsT_[(accT_>np.mean(accT_))*(accT_>=accA_)*(accT_>=accC_)*(accT_>=accG_)*(pred_==dDNA[1])] = 1
tfbsC_[(accC_>np.mean(accC_))*(accC_>=accA_)*(accC_>=accT_)*(accC_>=accG_)*(pred_==dDNA[2])] = 1
tfbsG_[(accG_>np.mean(accG_))*(accG_>=accA_)*(accG_>=accC_)*(accG_>=accT_)*(pred_==dDNA[3])] = 1
tfbs_ = np.array([tfbsA_,tfbsC_,tfbsT_,tfbsG_]) 
  
for i in range(1, len(pred_)):
    meanACCb[i] = np.mean(sum(pred_b[0:i] ==x[0:i])/(i+1))
    stdACCb[i] = np.std(meanACCb[0:i])
    meanACCf[i] = np.mean(sum(pred_f[0:i] ==x[0:i])/(i+1))
    stdACCf[i] = np.std(meanACCf[0:i])
    meanACC[i] = np.mean(sum(pred_[0:i] ==x[0:i])/(i+1))
    stdACC[i] = np.std(meanACC[0:i])
    
     
plt.figure
plt.subplot(411)
plt.plot(acc_.T)
plt.title('acc with mix L = 1 and L = 5')
plt.legend(['A','T','C','G'])
plt.subplot(412)
plt.plot(TP_.T)
plt.title('sensitivity with mix L = 1 and L = 5')
plt.legend(['A','T','C','G'])
plt.subplot(413)
plt.plot(TN_.T)
plt.title('specificity with mix L = 1 and L = 5')
plt.legend(['A','T','C','G'])
plt.subplot(414)
plt.plot(tfbs_.T)
plt.title('accordant tfbs with mix L = 1 and L = 5')
plt.legend(['A','T','C','G'])


plt.figure
plt.subplot(311)
plt.plot(range(len(meanACCb)),meanACCb,'k') 
plt.fill_between(range(len(meanACCb)),meanACCb-stdACCb,meanACCb+stdACCb)
plt.title('total ACC with L = 5')

plt.subplot(312)
plt.plot(range(len(meanACCf)),meanACCf,'k') 
plt.fill_between(range(len(meanACCf)),meanACCf-stdACCf,meanACCf+stdACCf)
plt.title('total ACC with L = 1')

plt.subplot(313)
plt.plot(range(len(meanACC)),meanACC,'k') 
plt.fill_between(range(len(meanACC)),meanACC-stdACC,meanACC+stdACC)
plt.title('total ACC with mix L = 1 and L = 5') 

         
#if max(0, k -1) == max(0, k - Lf):
#                nfn[j] = Seq[max(0, k-1)].count(LDNA[j]) + nfn[j]
#                nfn_[j] = len(x) - nfn[j] + nfn_[j]
#            else:
#                Min = min(max(0, k-1), max(0, k-Lf))
#                Max = max(max(0, k-1), max(0, k-Lf))
#                nfn[j] = Seq[Min:Max].count(LDNA(j)) + nfn[j]
#                nfn_[j] = len(x) - nfn[j] + nfn_[j]
#            if max(0, k) == max(0, k - Lf):
#                nfd[j] = Seq[max(0, k)].count(LDNA[j]) + nfd[j]
#                nfd_[j] = len(x)-nfd[j] + nfd_[j]
#            else:
#                Min = min(max(0, k), max(0, k-Lf))
#                Max = max(max(0, k), max(0, k-Lf))
#                nfd[j] = Seq[Min:Max].count(LDNA(j)) + nfd[j]
#                nfd_[j] = len(x)-nfd[j] + nfd_[j]
#            if max(0, k -1) == max(0, k - Lb):
#                nbn[j] = Seq[max(0, k-1)].count(LDNA[j]) + nbn[j]
#                nbn_[j] = len(x)-nbn[j] + nbn_[j]
#            else:
#                Min = min(max(0, k-1), max(0, k-Lb))
#                Max = max(max(0, k-1), max(0, k-Lb))
#                nbn[j] = Seq[Min:Max].count(LDNA(j)) + nbn[j]
#                nbn_[j] = len(x)-nbn[j] + nbn_[j]
#            if max(0, k) == max(0, k - Lb):
#                nbd[j] = Seq[max(0, k)].count(LDNA[j]) + nbd[j]
#                nbd_[j] = len(x)-nbd[j] + nbd_[j]
#            else:
#                Min = min(max(0, k), max(0, k-Lb))
#                Max = max(max(0, k), max(0, k-Lb))
#                nbd[j] = Seq[Min :Max].count(LDNA(j)) + nbd[j] 
#                nbd_[j] = len(x)-nbd[j] + nbd_[j]            
           
#
#    Seq[i] = max(s)
#    for j in range(len(LDNA)):
#        if max(0, i - 1) == max(0, i - Lb):
#            nbn[j] = Seq[max(0, i-1)].count(LDNA[j])
#            if nbn[j] == 0:
#                nbn[j] = Seq[max(0, i-1)].count(LDNA[j])
#            nbn_[j] = 1-nbn[j]
#        else:
#            Min = min(max(0, i-1), max(0, i-Lb))
#            Max = max(max(0, i-1), max(0, i-Lb))
#            nbn[j] = Seq[Min:Max].count(LDNA(j))
#            nbn_[j] = Lb-1-nbn[j]
#        if max(0, i) == max(0, i - Lb):
#            nbd[j] = Seq[max(0, i)].count(LDNA[j])
#            nbd_[j] = 1-nbd[j]
#        else:
#            Min = min(max(0, i), max(0, i-Lb))
#            Max = max(max(0, i), max(0, i-Lb))
#            nbd[j] = Seq[Min :Max].count(LDNA(j))
#            nbd_[j] = Lb-nbd[j]             
#       
#        Pf[j] = nfn[j]/(nfd[j]+0.0001)
#        Pf_[j] = nfn_[j]/(nfd_[j]+0.0001)
#            
#        if Pf[j]> Pf_[j]:
#            Lp0[j] = 1          
#        if  x[i] == LDNA[0]:
#            PF[0,j] = max(Pf[j], Pf_[j])
##            PB[0,j] = max(Pb[j], Pb_[j])               
#        if  x[i] == LDNA[1]: 
#            PF[1,j] = max(Pf[j], Pf_[j])
##            PB[1,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[2]: 
#            PF[2,j] = max(Pf[j], Pf_[j])
##            PB[2,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[3]: 
#            PF[3,j] = max(Pf[j], Pf_[j])
##            PB[3,j] = max(Pb[j], Pb_[j])     
#        if j == 3:
#            nfn_ = np.array(list([0]*len(LDNA)))
#            nfd_ = np.array(list([0]*len(LDNA)))
#            nfn = np.array(list([0]*len(LDNA)))
#            nfd = np.array(list([0]*len(LDNA)))
##            nbn_ = np.array(list([0]*len(LDNA)))
##            nbd_ = np.array(list([0]*len(LDNA)))
##            nbn = np.array(list([0]*len(LDNA)))
##            nbd = np.array(list([0]*len(LDNA))) 
#            Lp0 = np.array(list([0]*len(LDNA)))    
#            
#            
#        Pb[j] = nbn[j]/(nbd[j]+0.0001)
#        Pb_[j] = nbn_[j]/(nbd_[j]+0.0001)
#            
#        if Pf[j]> Pf_[j] & k == i:
#            Lp0[j] = 1          
#        if  x[i] == LDNA[0]:
#            PF[0,j] = max(Pf[j], Pf_[j])
#            PB[0,j] = max(Pb[j], Pb_[j])               
#        if  x[i] == LDNA[1]: 
#            PF[1,j] = max(Pf[j], Pf_[j])
#            PB[1,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[2]: 
#            PF[2,j] = max(Pf[j], Pf_[j])
#            PB[2,j] = max(Pb[j], Pb_[j])
#        if  x[i] == LDNA[3]: 
#            PF[3,j] = max(Pf[j], Pf_[j])
#            PB[3,j] = max(Pb[j], Pb_[j])  
#    
#        if j == 3:
#            nfn_ = np.array(list([0]*len(LDNA)))
#            nfd_ = np.array(list([0]*len(LDNA)))
#            nfn = np.array(list([0]*len(LDNA)))
#            nfd = np.array(list([0]*len(LDNA)))
#            nbn_ = np.array(list([0]*len(LDNA)))
#            nbd_ = np.array(list([0]*len(LDNA)))
#            nbn = np.array(list([0]*len(LDNA)))
#            nbd = np.array(list([0]*len(LDNA))) 
#            Lp0 = np.array(list([0]*len(LDNA)))
        
        
        
        
    
      
                