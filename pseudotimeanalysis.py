# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:05:13 2019

@author: LilyHeAsamiko
"""
 


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
seq3 = np.array(dataset)                                                           et4.iloc[:,0].astype('str'),dtype = 'str')
dataset4 = dataset4.T

# peeking at the dataset
dataset.head().T

dataset4 = dataset4.iloc[row,:]

stats1 = dataset1.describe()
dataset1 = dataset1.astype('float16')
tg_mean = np.mean(dataset1,0)
tg_std = np.std(dataset1,0)
dataset4 = dataset4.astype('float16')

def skewX(X):
    kg = []
    kmin = 0
    ti = 0
    if np.size(np.shape(X)) ==1:
        for i in range(1,np.size(X,0)):
            if i == 1:
                kmin = abs(X[i]-X[i-1])
            if X[i]-X[i-1] < kmin:
                kmin = abs(X[i]-X[i-1])
                ti = i-1
            kg.append(X[i]-X[i-1])
    return np.array(kg), kmin,ti

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

def grad(X, y):
    G = []

#    G =0
    for n in range(X.shape[1]):
        k,km,ti = skewX(X[:,n])
        print(['k:'+str(k)])
        tc = y
        tg = np.concatenate((np.array(range(tc)),np.array(range(tc+1,np.size(X,0)))),axis= 0)
        mu = np.mean(X[:,n])
        numerator = 2 * mu
        denominator = 1 + np.exp(-np.dot(k.T, tc-tg))+0.0001
        print(numerator)
        print(['d:'+str(denominator)])
        G.append(numerator / denominator)

    return G,km,ti

def grad1(w1_val, w2_val, y):
    G = 0

    for n in range(X.shape[1]):
        numerator = np.exp(-y[n] * np.dot(w, X[:,n])) * (-y[n]) * X[:,n]
        denominator = 1 + np.exp(-y[n] * np.dot(w, X[:,n]))

        G += numerator / denominator

    return G

def log_loss(k, X, y):
    L = []
    
#    for n in range(X.shape[1]):

    tc = y
#        tg = np.concatenate((np.array(range(tc)),np.array(range(tc+1,np.size(X,0)))),axis= 0)
    tg = range(np.size(X,1)) 
#    if len(X.shape) == 2:
#        Y = np.tile(y, (X.shape[0], 1)).T
#    else:
#        Y = y
    L.append(np.std(np.log(1 + np.exp(-k.T*(tc-tg)))))
#    L.append(np.std(np.log(1 + np.exp(-k*(tc-tg)))))

    L = np.mean(L, axis = 0)
    return L

def log_loss1(w1_val, w2_val, y):

    if len(w.shape) == len(y):
        Y = np.tile(y, (w.T.shape[1], 1)).T
    else:
        Y = y

    const = -Y * np.dot(X, w.T)
    L = np.sum(np.log(1 + np.exp(const)), axis = 0)
    return L

def regression(w0_val, w0, w1_val, w1, dataset4, seq3,  cell_name, tc, stepsize, steps, i):
    # latent t = 0, test the distribution 
    cols_0, cols_1 = test_normal(w0_val, seq3)
    
    # cell lines
    tc_lb, tc_err = cell2digit(tc[:,1], cell_name)
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
    T = np.size(w0_val,0)
    R = np.size(w2_val,0)
    C = np.size(w2_val,1) 
    X = np.zeros((T,R,C))
    # normalize time-series  ///////////   
    for c in range(C):
        w0_val[:, c] = (w0_val[:, c] - w0_val[:, c].mean())/w0_val[:, c].std()

    for i in range(T):
        X[i,:,:] = w2_val*np.repeat(w0_val[i,:],R,axis = 0).reshape(R,C)
 
    y, y_err = cell2digit(np.array(cell[cell_name[tc_lb_1]].index, dtype = 'str'), cell, unique = 1)
    y = np.array(y)
       
#    y = np.array(tc[0:9,2].astype('float32'))
#    w2,w2_val = data2cell(w1, tc_lb_1)
    #X = np.dot(w1_val,w2_val.T)


    ## latent t = tl
#    tl = 5
  
#        # Predict cell line A1 probability      
    for r in range(R): 
    #test the regressed data for cell_r
        Cols_0 = []
        Cols_1 = []
        Cols_0[r,:], Cols_1[r,:] = test_normal(Kcell[:,r,:], seq3)                         
        tc = y[r]
#       tg = np.concatenate((np.array(range(tc)),np.array(range(tc+1,np.size(X,0)))),axis= 0)
        tg = np.array(range(T))           
        y_1 = tc[range(R+1),2]
        y_b = (y_1 > 0).astype(int)
        y_b = 2*y_b - 1
            
        k_pred = np.linspace(-1, 1, C)   
        K_pred = []                                                                    
    #    
    #    X = X[:,idx]
    #    y = y[idx]
        #    W = []
        accuracies = []
        losses_pred = []
#            E = []            
        accuracy = 0
        P = []
        p = []
        for iteration in range(steps):
            dk_pred,km_pred,ti = grad(Kcell[:,r,:]*1000, tc)
            tc = ti
            step_size = Kcell[tc,r,:]
            dk_pred,km_pred,ti = grad(Kcell[range(tc),r,:]*1000, tc)
            k_pred = k_pred - step_size * np.array(dk_pred)
#                k = k - step_size * dk
            K_pred.append(k_pred)
            loss_val = log_loss(k_pred, Kcell[range(tc),r,:]*1000, tc)
            print (loss_val)
            print ("Iteration %d: k = %s (log-loss = %.2f)" % \
                  (iteration, str(k_pred), loss_val))

            losses_pred.append(loss_val)
            mu = Kcell[tc,r,:].mean()
            y_prob = 2*mu / (1 + np.exp(np.dot(np.dot(-k_pred.T,*Kcell[:,r,:].T),(tc-tg))))
#    #           # Threshold at 0.5 (results are 0 and 1)
            y_pred = (y_prob > 0.5).astype(int)
#                    # Transform [0,1] coding to [-1,1] coding
            y_pred = 2*y_pred - 1

            accuracy += y_pred == y_b[r]
            accuracies.append(accuracy/iteration)

            if iteration == 0:
                l0 = loss_val
                lmin = l0
                iterminl = 0
            if loss_val < l0:
                lmin = loss_val
                iterminl = iteration
                p = 2*mu / (1 + np.exp(np.dot(-k_pred.T,*Kcell[:,r,:].T).T*(tc-tg)))
            
        losses_pred = np.array(losses_pred)
        accuracies =np.array(accuracies)
        P.append(np.array(p))
            
        plt.figure(r)
        plt.subplot(2,1,1)
        plt.plot(range(steps), losses_pred[range(steps)])        
        plt.title(['A1_'+str(cell_1[r])+'prediction of tc with least average mse loss at Best epoch'])
        plt.subplot(2,1,2)
        plt.plot(range(steps), accuracies[range(steps)])
        plt.title(['A1_'+str(cell_1[r])+'prediction accuracy at Best epoch']) 
        plt.tight_layout
        plt.savefig(str(['A1_'+str(cell_1[r])+'_50_prediction']),format= 'png', pdi = 200)
            
        plt.figure(r)
        plt.subplot(3,1,1)
        plt.scatter(range(C), Kcell[tc,r,:]*1000)        
        plt.title(['A1_'+str(cell_1[r])+'expression of pseudotime'])
        plt.subplot(3,1,2)
        plt.scatter(range(T), Kcell[range(T),r,:]*1000)
        plt.title(['A1_'+str(cell_1[r])+'expression evolution']) 
        plt.subplot(3,1,3)
        plt.scatter(range(T), p)
        plt.title(['A1_'+str(cell_1[r])+'cell evolution'])             
        plt.tight_layout
        plt.savefig(str(['A1_'+str(cell_1[r])+'_50_expression']),format= 'png', pdi = 200)
            
            
#    #
#    W = np.array(W)
#    Y = np.dot(W,X)
    #    emin = np.min(E)
    #    # Plot the path.
    #
    #    xmin = min(np.min(X,axis = 0))
    #    xmax = max(np.max(X,axis = 0))
    #    ymin = min(np.min(X,axis = 1))
    #    ymax = max(np.max(X,axis = 1))
    #
    #    Xg, Yg = np.meshgrid(np.linspace(xmin, xmax, np.size(W,0)),
    #                         np.linspace(ymin, ymax, np.size(W,1)))
    ##    Zg = np.reshape(Zg, Xg.shape)
    ##    levels = np.linspace(70, 100, 20)
    #
    #    fig = plt.figure()
    ##    ax = fig.add_subplot(111, projection='3d')    
    #    plt.pcolor(Xg, Yg, W.T,
    #                  alpha=0.5,
    #                  cmap=plt.cm.bone)
    #
    #    fig = plt.figure()
    #    plt.plot(W[itermin,:], 'ro-')
    #    plt.xlabel('w$_0$')
    #    plt.ylabel('w$_1$')
    #    plt.title('Optimization path')
    #    plt.grid()
    #    plt.axis([xmin, xmax, ymin, ymax])
    #
    ##    plt.annotate('Starting point',
    ##             xy=(W[0, 0], W[0, 1]),
    ##             xytext=(-0.5, 0),
    ##             size=13,
    ##             bbox=dict(boxstyle="round4", fc="w", ec = "g"),
    ##             arrowprops=dict(arrowstyle="simple",
    ##                             connectionstyle="arc3,rad=0.15",
    ##                             shrinkA = 0,
    ##                             shrinkB = 8,
    ##                             fc = "g",
    ##                             ec = "g"))
    ##
    ##    plt.annotate('Endpoint',
    ##             xy=(W[-1, 0], W[-1, 1]),
    ##             xytext=(1.2, 0),
    ##             size=13,
    ##             bbox=dict(boxstyle="round4", fc="w", ec = "g"),
    ##             arrowprops=dict(arrowstyle="simple",
    ##                             connectionstyle="arc3,rad=0.15",
    ##                             shrinkA = 0,
    ##                             shrinkB = 8,
    ##                             fc = "g",
    ##                             ec = "g"))    
    #    
    #    plt.show()
    #             
    #    plt.figure() 
    #    plt.plot(100.0 * np.array(accuracies), linewidth = 2,
    #               label = "Classification Accuracy")
    #    plt.ylabel('Accuracy / %')
    #    plt.xlabel('Iteration')
    #    plt.legend(loc = 4)
    #    plt.grid()
    #    plt.annotate('best_epoch', xy = (itermin,accuracies[itermin])) 
    #    plt.show()
    #    plt.tight_layout()
    #    plt.savefig(['log_loss_minimization'+str(i)+'.pdf'], bbox_inches = "tight")
    #    
    #    plt.figure()        
    #    plt.plot(steps, accuracies)
    #    plt.savefig(['accuracy'+str(i)+'.pdf'], bbox_inches = "tight")
    #
    #    
    #    return W[itermin,:],Y[itermin,:],accuracies,emin

#    for col in idx:
#        w1_val[:, col] = w1_val[:, col] - w1_val[:, col].mean()
#    
#    X = X[:,idx]
#    y = y[idx]
#
#    step_size = stepsize
#    W = []
#    accuracies = []
#    losses = []
#    E = []
#    
#    for iteration in range(steps):
#        w = w - step_size * grad(w, X, y)
#
#        loss_val = log_loss(w, X, y)
#        print (loss_val)
#        print ("Iteration %d: w = %s (log-loss = %.2f)" % \
#              (iteration, str(w), loss_val))
#
#        # Predict class 1 probability
#        y_prob = 1 / (1 + np.exp(-np.dot(w.T, X)))
#                # Threshold at 0.5 (results are 0 and 1)
#        y_pred = (y_prob > 0.5).astype(int)
#                # Transform [0,1] coding to [-1,1] coding
#        y_pred = 2*y_pred - 1
#        y_b = (y > 0).astype(int)
#        y_b = 2*y_b - 1
#        
#        accuracy = np.mean(y_pred == y_b)
#        accuracies.append(accuracy)
#        losses.append(loss_val)
#        if iteration == 0:
#            e0 = y-w*X
#            emin = e0
#            itermin = 0
#        e = y-w*X
#        if e < e0:
#            emin = e
#            itermin = iteration
#        E.append(e)
#        W.append(w)
#
#    W = np.array(W)
#    Y = np.dot(W,X)
##    emin = np.min(E)
#    # Plot the path.
#
#    xmin = min(np.min(X,axis = 0))
#    xmax = max(np.max(X,axis = 0))
#    ymin = min(np.min(X,axis = 1))
#    ymax = max(np.max(X,axis = 1))
#
#    Xg, Yg = np.meshgrid(np.linspace(xmin, xmax, np.size(W,0)),
#                         np.linspace(ymin, ymax, np.size(W,1)))
##    Zg = np.reshape(Zg, Xg.shape)
##    levels = np.linspace(70, 100, 20)
#
#    fig = plt.figure()
##    ax = fig.add_subplot(111, projection='3d')    
#    plt.pcolor(Xg, Yg, W.T,
#                  alpha=0.5,
#                  cmap=plt.cm.bone)
#
#    fig = plt.figure()
#    plt.plot(W[itermin,:], 'ro-')
#    plt.xlabel('w$_0$')
#    plt.ylabel('w$_1$')
#    plt.title('Optimization path')
#    plt.grid()
#    plt.axis([xmin, xmax, ymin, ymax])
#
##    plt.annotate('Starting point',
##             xy=(W[0, 0], W[0, 1]),
##             xytext=(-0.5, 0),
##             size=13,
##             bbox=dict(boxstyle="round4", fc="w", ec = "g"),
##             arrowprops=dict(arrowstyle="simple",
##                             connectionstyle="arc3,rad=0.15",
##                             shrinkA = 0,
##                             shrinkB = 8,
##                             fc = "g",
##                             ec = "g"))
##
##    plt.annotate('Endpoint',
##             xy=(W[-1, 0], W[-1, 1]),
##             xytext=(1.2, 0),
##             size=13,
##             bbox=dict(boxstyle="round4", fc="w", ec = "g"),
##             arrowprops=dict(arrowstyle="simple",
##                             connectionstyle="arc3,rad=0.15",
##                             shrinkA = 0,
##                             shrinkB = 8,
##                             fc = "g",
##                             ec = "g"))    
#    
#    plt.show()
#             
#    plt.figure() 
#    plt.plot(100.0 * np.array(accuracies), linewidth = 2,
#               label = "Classification Accuracy")
#    plt.ylabel('Accuracy / %')
#    plt.xlabel('Iteration')
#    plt.legend(loc = 4)
#    plt.grid()
#    plt.annotate('best_epoch', xy = (itermin,accuracies[itermin])) 
#    plt.show()
#    plt.tight_layout()
#    plt.savefig(['log_loss_minimization'+str(i)+'.pdf'], bbox_inches = "tight")
#    
#    plt.figure()        
#    plt.plot(steps, accuracies)
#    plt.savefig(['accuracy'+str(i)+'.pdf'], bbox_inches = "tight")
#
#    
#    return W[itermin,:],Y[itermin,:],accuracies,emin

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
ax.bar(xx,C,yy)
ax = fig.add_subplot(212, projection='3d')
ax.scatter(xx,yy,C)
plt.show()
plt.savefig('cell evolution_3D.pdf',dpi=200)

w1,w1_val = data2cell(dataset1, lb)
w0,w0_val = data2cell(dataset4, lb)




#
##standardize the data to normal distribution
#dataset1_std = preprocessing.scale(dataset1)
#dt1_df = pd.DataFrame(dataset1_std.astype(np.float16))
##replace blank cells with mean
#dt1_df.fillna(np.mean(np.mean(dt1_df)), inplace=True)
#dt1_df.dropna(how = 'all')
#dataset1 = dt1_df.values

