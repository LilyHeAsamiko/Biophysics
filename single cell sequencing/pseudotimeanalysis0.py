# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 21:29:55 2019

@author: LilyHeAsamiko
"""
import pandas as pd
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import unicodedata
import sys
sys.path.append('D:/TUT/Medical/biophysics/experiment')
import math
test = math.inf


def data2cell(dataset, lb):
    rows = lb
    w1_val = []
    for r in rows:
#        w1.append(dataset.iloc[:,c+1])
        w1_val.append(dataset.values[r,:].astype('float32'))
    return dataset.iloc[rows,:],np.array(w1_val)

def cell2digit(col, cell_name, unique = 0):
#    col_N = np.unique(col)
    N = len(cell_name)
    label = []
    notfound = []
    ind = range(N)
    i = 0
    if unique == 1:
        cell_name1, ind = np.unique(cell_name, return_index = True)       
        N = len(cell_name1)
    for c in col: 
        for n in range(N):
            if unique == 0:
                if cell_name[n].casefold() == c.casefold():
                    label.append(ind[n])
                elif n == ind[-1]:
                    notfound.append([i,c])
                    print("{:s} ".format(c)+ 'cannot be found labeled at'+str(i))
            else:
                if cell_name1[n].casefold() == c.casefold():
                    label.append(ind[n])
                elif n == ind[-1]:
                    notfound.append([i,c])
                    print("{:s} ".format(c)+ 'cannot be found labeled at'+str(i))                
        i += 1
    return label,notfound

def name2digit(col):
    col_N = np.unique(col)
    N = len(col_N)
    label = []
    for c in col:
        for n in range(0, N):            
            if c.casefold() == col_N[n].casefold():
                label.append(n)
    return label

def test_normal(data, seq):
    cols = np.size(data,1)
    cols_0 = []
    cols_1 = []
    for c in range(0, cols):   
        ku0, sk0 = stats.normaltest(data[:,c])
        alpha = 0.05
        if sk0 < alpha:  # null hypothesis: x comes from a normal distribution
            print("for gene {:s}, p = {:g}".format(seq[c], sk0)+'The null hypothesis: normal distribution is not followed can be rejected')
            cols_1.append(c)
        else:
            print('\x1b[1;34;43m'+"for gene: {:s} at t0:{:g} , p = {:g} ".format(seq[c], c, sk0)+"The null hypothesis cannot be rejected"+ '\x1b[0m')
            cols_0.append(c)    
    return cols_0, cols_1 

def grad(w, X, y):

    G = 0

    for n in range(X.shape[0]):
        if np.size(X.shape) == 2:
            numerator = np.exp(-y[n] * np.dot(w, X[n])) * (-y[n]) * X[n]
            denominator = 1 + np.exp(-y[n] * np.dot(w, X[n]))
        else:
            numerator = np.exp(-y[n] * np.dot(w, X[n])) * (-y[n]) * X[n]
            print(['n:'+ str(numerator)])
            denominator = 1 + np.exp(-y[n] * np.dot(w, X[n]))
            print(['d:'+ str(denominator)])
    
        G += numerator / denominator

    return G

def log_loss(w, X, y):

    if len(w.shape) == 2:
        Y = np.tile(y, (w.T.shape[1], 1)).T
    else:
        Y = y
#    print(['product: ', str(np.dot(X,w))])
#    print(['exp: ', str(-Y * np.dot(X,w))])
#    print(['log: ', str(np.log(1 + np.exp(-Y * np.dot(X,w))))])
    L = np.mean(np.log(1 + np.exp(-Y * np.dot(X,w))), axis = 0)
    return L

def loss (w, X, y):

    if len(w.shape) == 2:
        Y = np.tile(y, (w.T.shape[1], 1)).T
    else:
        Y = y

    L = np.mean(1 + np.exp(-Y * np.dot(X,w)), axis = 0)
    return L
#def regression(w0_val, w0, w1_val, w1, dataset4, seq3,  cell_name, Mc, stepsize, steps, i):
#        
##    y = np.array(tc[0:9,2].astype('float32'))
##    w2,w2_val = data2cell(w1, tc_lb_1)
#    #X = np.dot(w1_val,w2_val.T)
#
#
#    ## latent t = tl
##    tl = 5
#    T_ind = reversed(range(T))
#    Lcell = np.zeros((T, R, C))
#    Kcell = np.zeros((T, R, C))
#    
#    for t in T_ind:
#        # regression on different cell_r of interests
#        for r in range(R):                        
#            tc = y[r]
##            tg = np.concatenate((np.array(range(tc)),np.array(range(tc+1,np.size(X,0)))),axis= 0)
#            tg = np.array(range(t+1))
#    
#            k = np.linspace(-1, 1, C)   
#            K = []                                                                    
#    #    
#    #    X = X[:,idx]
#    #    y = y[idx]
#            step_size = w0_val[t,:]
#        #    W = []
##            accuracies = []
#            losses = []
#            E = []            
#        #    
#            for iteration in range(steps):
#                dk,km,ti = grad(X[range(t+1),r,:]*1000, tc)
#                k = k - step_size * np.array(dk)
##                k = k - step_size * dk
#                K.append(k)
#                loss_val = log_loss(k, X[range(t+1),r,:]*1000, tc)
#                print (loss_val)
#                print ("Iteration %d: k = %s (log-loss = %.2f)" % \
#                       (iteration, str(k), loss_val))
#
#                losses.append(loss_val)
#                e = np.std(tc-np.dot(k.T, w0_val[range(t+1)].T)*tg.T)
##                e = tc-k*t
#                E.append(e)
#                if iteration == 0:
#                    e0 = np.std(tc-np.dot(k.T, w0_val[range(t+1)].T)*tg.T)
##                    e0 = tc-k*t
#                    emin = e0
#                    itermine = 0
#                    l0 = loss_val
#                    lmin = l0
#                    iterminl = 0
#                if loss_val < l0:
#                    lmin = loss_val
#                    iterminl = iteration
#                if e < e0:
#                    emin = e
#                    itermine = iteration
#                    
#            E = np.array(E)
#            losses = np.array(losses)
#            K = np.array(K)
#            plt.figure(r)
#            plt.subplot(2,1,1)
#            plt.plot(range(steps), losses[range(steps)])        
#            plt.title(['A1_'+str(cell_1[r])+'least average mse loss at Best epoch'])
#            plt.subplot(2,1,2)
#            plt.plot(range(steps), E[range(steps)])
#            plt.title(['A1_'+str(cell_1[r])+'least average mse residue at Best epoch']) 
#            plt.tight_layout
#            plt.savefig(str(['A1_'+str(cell_1[r])+'_50']),format= 'png', pdi = 200)
#            
#            Kcell[t,r,:] = np.array(K[iterminl,])
#            Lcell[t,r,:] = np.array(losses[iterminl,])
#
##        # Predict cell line A1 probability      
#    for r in range(R): 
#    #test the regressed data for cell_r
#        Cols_0 = []
#        Cols_1 = []
#        Cols_0[r,:], Cols_1[r,:] = test_normal(Kcell[:,r,:], seq3)                         
#        tc = y[r]
##       tg = np.concatenate((np.array(range(tc)),np.array(range(tc+1,np.size(X,0)))),axis= 0)
#        tg = np.array(range(T))           
#        y_1 = tc[range(R+1),2]
#        y_b = (y_1 > 0).astype(int)
#        y_b = 2*y_b - 1
#            
#        k_pred = np.linspace(-1, 1, C)   
#        K_pred = []                                                                    
#    #    
#    #    X = X[:,idx]
#    #    y = y[idx]
#        #    W = []
#        accuracies = []
#        losses_pred = []
##            E = []            
#        accuracy = 0
#        P = []
#        p = []
#        for iteration in range(steps):
#            dk_pred,km_pred,ti = grad(Kcell[:,r,:]*1000, tc)
#            tc = ti
#            step_size = Kcell[tc,r,:]
#            dk_pred,km_pred,ti = grad(Kcell[range(tc),r,:]*1000, tc)
#            k_pred = k_pred - step_size * np.array(dk_pred)
##                k = k - step_size * dk
#            K_pred.append(k_pred)
#            loss_val = log_loss(k_pred, Kcell[range(tc),r,:]*1000, tc)
#            print (loss_val)
#            print ("Iteration %d: k = %s (log-loss = %.2f)" % \
#                  (iteration, str(k_pred), loss_val))
#
#            losses_pred.append(loss_val)
#            mu = Kcell[tc,r,:].mean()
#            y_prob = 2*mu / (1 + np.exp(np.dot(np.dot(-k_pred.T,*Kcell[:,r,:].T),(tc-tg))))
##    #           # Threshold at 0.5 (results are 0 and 1)
#            y_pred = (y_prob > 0.5).astype(int)
##                    # Transform [0,1] coding to [-1,1] coding
#            y_pred = 2*y_pred - 1
#
#            accuracy += y_pred == y_b[r]
#            accuracies.append(accuracy/iteration)
#
#            if iteration == 0:
#                l0 = loss_val
#                lmin = l0
#                iterminl = 0
#            if loss_val < l0:
#                lmin = loss_val
#                iterminl = iteration
#                p = 2*mu / (1 + np.exp(np.dot(-k_pred.T,*Kcell[:,r,:].T).T*(tc-tg)))
#            
#        losses_pred = np.array(losses_pred)
#        accuracies =np.array(accuracies)
#        P.append(np.array(p))
#            
#        plt.figure(r)
#        plt.subplot(2,1,1)
#        plt.plot(range(steps), losses_pred[range(steps)])        
#        plt.title(['A1_'+str(cell_1[r])+'prediction of tc with least average mse loss at Best epoch'])
#        plt.subplot(2,1,2)
#        plt.plot(range(steps), accuracies[range(steps)])
#        plt.title(['A1_'+str(cell_1[r])+'prediction accuracy at Best epoch']) 
#        plt.tight_layout
#        plt.savefig(str(['A1_'+str(cell_1[r])+'_50_prediction']),format= 'png', pdi = 200)
#            
#        plt.figure(r)
#        plt.subplot(3,1,1)
#        plt.scatter(range(C), Kcell[tc,r,:]*1000)        
#        plt.title(['A1_'+str(cell_1[r])+'expression of pseudotime'])
#        plt.subplot(3,1,2)
#        plt.scatter(range(T), Kcell[range(T),r,:]*1000)
#        plt.title(['A1_'+str(cell_1[r])+'expression evolution']) 
#        plt.subplot(3,1,3)
#        plt.scatter(range(T), p)
#        plt.title(['A1_'+str(cell_1[r])+'cell evolution'])             
#        plt.tight_layout
#        plt.savefig(str(['A1_'+str(cell_1[r])+'_50_expression']),format= 'png', pdi = 200)

if __name__ == "__main__":

plt.close("all")

path = r'D:\TUT\Medical\biophysics\experiment'
#importing the dataset
#dataset2 = np.loadtxt('D:/TUT/Medical/biophysics/experiment/TF predictivity.txt', delimiter ='\t',dtype = 'str')
dataset = pd.read_csv(f'{path}/TF predictivity.txt',sep='\t')
seq1 = np.array(dataset.iloc[:,0].astype('str'),dtype = 'str')
#dataset2 = np.loadtxt('D:/TUT/Medical/biophysics/experiment/microarray sources.txt', delimiter ='\t',dtype = 'str')
dataset2 = pd.read_csv(f'{path}/Partially reprogrammed cells Z-score.txt',sep='\t')
seq2 = np.array(dataset2.iloc[:,0].astype('str'),dtype = 'str')
col = range(0,3)
row = range(1,np.size(dataset,1))
dataset3 = pd.read_csv(f'{path}/microarray sources.txt',sep='\t', usecols = col)
dataset1 = dataset.T.iloc[row,:]
cell_name = np.array(dataset1.index.astype('str'),dtype = 'str')
dataset4 = pd.read_csv(f'{path}/TF Z-score.txt',sep='\t')
seq3 = np.array(dataset4.iloc[:,0].astype('str'),dtype = 'str')
dataset4 = dataset4.T

# peeking at the dataset
dataset.head().T

dataset4 = dataset4.iloc[row,:]

stats1 = dataset1.describe()
dataset1 = dataset1.astype('float16')
tg_mean = np.mean(dataset1,0)
tg_std = np.std(dataset1,0)
dataset4 = dataset4.astype('float16')

GSE = dataset3.iloc[:,0].astype('str')
GSE_stat = GSE.describe()
GSM = dataset3.iloc[:,1].astype('str')
GSM_stat = GSM.describe() 
cell = dataset3.iloc[:,2].astype('str')
cell_stat = cell.describe() 

Mc = np.zeros((51,3))
Mc= np.array([['A2','ESC',0.178],['A2','MSC',0.158],['A2','myoblast',0.142],['A2','Megakaryocyte-Erythroid Progenitor (MEP)',0.129],['A2','endothelial cell - blood vessel',0.113],['A2','keratinocyte',0.112],['A2','medullary thymic epithelial',-0.111],['A2','adipose - brown',-0.117],['A2','natural killer cells',-0.13],['A2','Common Myeloid Progenitor (CMP)',-0.138],
['B3','ESC',0.222],['B3','MSC',0.161],['B3','endothelial cell - blood vessel', 0.139],['B3','myoblast',0.138],['B3','Granulocyte-Monocyte Progenitor (GMP)',0.127],['B3','kidney',0.111],['B3','Megakaryocyte-Erythroid Progenitor (MEP)',-0.107],['B3','cornea',0.107],['B3','natural killer cells',-0.129],
['B1V1+','myoblast',0.181],['B1V1+','prostate',0.164],['B1V1+','MSC',0.154],['B1V1+','Megakaryocyte-Erythroid Progenitor (MEP)',0.138],['B1V1+','keratinocyte',0.136],['B1V1+','cornea',0.125],['B1V1+','ESC',0.111],['B1V1+','intestine - paneth cell',-0.111],['B1V1+','Common Myeloid Progenitor (CMP)',-0.122],
['B1V1-','ESC',0.382],['B1V1-','EpiSC',0.184],['B1V1-','Megakaryocyte-Erythroid Progenitor (MEP)',0.160],['B1V1-','myoblast',0.145],['B1V1-','NSC',-0.108],['B1V1-','T cell',0.115],['B1V1-','skeletal muscle',-0.117],['B1V1-','Common Myeloid Progenitor (CMP)',-0.154],
['MCV6','Megakaryocyte-Erythroid Progenitor (MEP)',0.155],['MCV6','myoblast',0.150],['MCV6','ESC',0.149],['MCV6','keratinocyte',0.145],['MCV6','Common Lymphoid Progenitor (CLP)',0.107],['MCV6','Granulocyte-Monocyte Progenitor (GMP)',0.107],['MCV6','cornea',0.107],['MCV6','Common Myeloid Progenitor (CMP)',-0.130],
['MCV8','ESC',0.203],['MCV8','Megakaryocyte-Erythroid Progenitor (MEP)',0.191],['MCV8','myoblast',0.160],['MCV8','cornea',0.119],['MCV8','prostate',0.113],['MCV8','skeletal muscle', -0.141],['MCV8','Common Myeloid Progenitor (CMP)',-0.142]])

lb, lb_err = cell2digit(np.array(cell, dtype = 'str'), cell_name)
lb = np.array(lb)
lb_b = preprocessing.Binarizer(lb)
xx = np.array(name2digit(GSE))
yy = np.array(name2digit(GSM))
C = lb/max(lb)
N = len(C)
[XX,YY] = np.sort(np.meshgrid(xx,yy))
CC = np.sort(np.meshgrid(C,C))
fig = plt.figure()
plt.scatter(xx,yy,CC)
plt.show()
plt.savefig('cell evolution.pdf',dpi=200)

fig = plt.figure()
ax = fig.add_subplot(211, projection='3d')
plt.title('cell_evolution: xaxis: cell, yaxis: expression, zaxis: psudotime')
ax.bar(xx,C,yy)
ax = fig.add_subplot(212, projection='3d')
ax.scatter(xx,yy,C)
plt.title('cell_evolution scatter on x_z_y')
plt.show()
plt.savefig('cell evolution_3D.pdf',dpi=200)

w1,w1_val = data2cell(dataset1, lb)
w0,w0_val = data2cell(dataset4, lb)

# latent t = 0, test the distribution 
cols_0, cols_1 = test_normal(w0_val, seq3)
    
# cell lines
tc_lb, tc_err = cell2digit(Mc[:,1], cell_name)
tc_lb_1 = np.array(tc_lb[0:9])
tc_lb_2 = np.array(tc_lb[10:18])
tc_lb_3 = np.array(tc_lb[19:27])
tc_lb_4 = np.array(tc_lb[28:35])
tc_lb_5 = np.array(tc_lb[36:43])
tc_lb_6 = np.array(tc_lb[44:50])


#    
#    idx = np.arange(y.size)
#    np.random.shuffle(idx)

#A1  
w2_val = np.array(dataset4.iloc[tc_lb_1,:].astype('float32'))
cell_1 = cell_name[tc_lb_1]
#    T = np.size(w0_val,0)
C = np.size(w2_val,0)
G = np.size(w2_val,1) 

# normalize time-series  
for g in range(G):
    w0_val[:, g] = (w0_val[:, g] - w0_val[:, g].mean())/w0_val[:, g].std()

#    for i in range(T):
#        X[i,:,:] = w2_val*np.repeat(w0_val[i,:],R,axis = 0).reshape(R,C)
#     
xc, xc_err = cell2digit(np.array(cell[cell_name[tc_lb_1]].index, dtype = 'str'), cell, unique = 1)
X = w0_val[xc,:]
y = np.array(Mc[tc_lb_1,2], dtype = float)
    

#    idg = np.zeros((1000,G))
#    for n in np.array(range(1000),dtype = int):
#        idg[n,:] = np.array(range(G),dtype = int).T
#        idg = np.array(range(G),dtype = int)
#        np.random.shuffle(idg[n,:])

#        X = X[:,idg[n,:]]
        
w = np.linspace(1,-1, G)
        
#step_size = np.mean(w2_val,axis =0)
step_size = 0.01
        
W = []
accuracies = []
losses = []
steps = 20
nmin = 0
idg = np.zeros((100,G))
for n in np.array(range(100),dtype = int):
    idg[n,:] = np.array(range(G),dtype = int).T
#    idg = np.array(range(G),dtype = int)
    np.random.shuffle(idg[n,])
    idg = np.array(idg, dtype = int)

    X = X[:,idg[n,:]]
    itermin = 0
    
    for iteration in range(steps):
    
        step_size = np.mean(w2_val,axis =0)
        w = w - step_size * grad(w, X, y)
        w = np.nan_to_num(w)     
        loss_val = log_loss(w, X, y)
        loss_val = np.nan_to_num(loss_val)  
        print (loss_val)
        print ("Iteration %d: w = %s (log-loss = %.2f)" % \
               (iteration, str(w), loss_val))
        
        # Predict class 1 probability
        y_prob = 1 / (1 + np.exp(-np.dot(X, w)))
        # Threshold at 0.5 (results are 0 and 1)
        y_pred = (y_prob > 0.5).astype(int)
        # Transform [0,1] coding to [-1,1] coding
        y_pred = 2*y_pred - 1
        # convert y
        y = (y > 0).astype(int)
        y = 2*y-1
        
        accuracy = np.mean(y_pred == y)
        accuracies.append(accuracy)
        losses.append(loss_val)
        
        W.append(w)
        
        if iteration == 0:
            l0 = loss_val
            lmin = l0

        if  loss_val < l0:
            lmin = loss_val
            itermin = iteration
            nmin = n
            wbest = w
            
        
W = np.array(W)

plt.figure,
plt.subplots_adjust(hspace=0.4, wspace=0.4)
plt.subplot(2,1,1)
plt.plot(range(len(W)))    
plt.title('losses')    
#plt.title(['A1_'+str(cell_1[r])+'least average mse loss at Best epoch'])
plt.subplot(2,1,2)
plt.plot(range(len(W)), accuracies)
plt.title('accuracy')    
#plt.title(['A1_'+str(cell_1[r])+'least average mse residue at Best epoch']) 
plt.tight_layout
#plt.savefig(str(['A1_'+str(cell_1[r])+'_50']),format= 'png', pdi = 200)
Xiter = X[:,idg[nmin]]
y_p_best = 1 / (1 + np.exp(-np.dot(Xiter, wbest)))


#            p = 2*mu / (1 + np.exp(np.dot(-k_pred.T,*Kcell[:,r,:].T).T*(tc-tg)))

#for i in np.array(range(C),dtype=int):
#    plt.figure,
#    plt.plot(np.array(range(G),dtype=int),Xiter[i,:])
#    plt.title('start TF')
acs = np.sort(accuracies)
idex, id_err = cell2digit(np.array(acs,dtype=str), np.array(accuracies,dtype=str))

Y = np.dot(W,X.T)
col = [idex[771466], idex[771469], idex[771560], idex[771564], idex[771565]]
Yb = Y[col]

for i in np.array(range(Yb.shape[0])):
    plt.figure(i),
    plt.title(['cell_best_predicted'])
    plt.scatter(range(Yb.shape[1]),Yb[i])
    plt.scatter(range(Yb.shape[1]),np.log(Yb[i]))
    plt.legend(['cell_best_predicted','log'])

#get pseudotime for cell
for i in np.array(range(C),dtype=int):
    fig = plt.figure(i)
    ax = fig.add_subplot(111, projection='3d')
    Cell_t = np.dot(w0_val, Xiter.T)
    Cell_T = Cell_t[:,i]
    ax.scatter(xx,yy,Cell_T)  
    plt.xlabel('pseudotime')
    plt.title(['expression for cell ',cell_1[i]])
    plt.show()
#check distribution
gene_pred = np.dot(Cell_t,Xiter);
test_normal(gene_pred.T, seq3)
A1_TF = np.dot(Cell_t.T,w0_val)
plt.figure(),
plt.subplot(211)
plt.imshow(Xiter,aspect = 'auto')
plt.title('TF_start')
plt.tight_layout()
plt.subplot(212)
plt.imshow(A1_TF,aspect = 'auto')
plt.title('TF_prediction_at cell line A1 ')
plt.tight_layout()