# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 09:38:21 2019

@author: user
"""

import csv
import sys
sys.path.append('D:/TUT/Medical/biophysics/experiment')
import os
if not os.path.exists('D:/TUT/Medical/biophysics/experiment'):
    os.makedirs('D:/TUT/Medical/biophysics/experiment')
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
#from keras.datasets import mnist
from keras import backend as K
from keras.models import Sequential, load_model
from keras.layers import Dense, Conv2D, MaxPool2D, Flatten, Reshape
from keras import optimizers
from keras.layers.core import Activation
from dpGANs import BasicDpGAN, DpGAN,  run
from sklearn.model_selection import train_test_split
import argparse

np.seterr(divide='ignore', invalid='ignore')
pm_paths = (['D:/TUT/Medical/biophysics/experiment/DpGAN_G.h5','D:/TUT/Medical/biophysics/experiment/DpGAN_D.h5'])
 

def parse_args():
    parser = argparse.ArgumentParser(prog='PROG', conflict_handler='resolve')
    parser.add_argument("-e", "--epochs", type=int, default=1000, help="the training epochs (set to 0 if you want to test dp-GAN only")     
    parser.add_argument("-se", "--sess", type=int, default=0, help="the training model (set to 0:{'GPU': 0}) If you wanna implement this code with GPU, change 0 to 1 or 2")
    parser.add_argument("-d", "--zdim", type=int, default=1000, help="the dmension of z")
    parser.add_argument("-is", "--img_size", type=int, default=26, help="the size of image")
    parser.add_argument("-l", "--lam", type=float, default=10., help="the improved WGAN regularization coefficient")
    parser.add_argument("-t", "--critic", type=int, default=4, help="the number of discriminator iterations per generator iteration")
    parser.add_argument("-m", "--batch", type=int, default=64, help="the batch size")
    parser.add_argument("-A", "--adam", type=float, nargs=3, default=(0.002, 0.5, 0.9), help="the adam optimizer parameters (alpha, beta1, beta2)")
    parser.add_argument("-C", "--clip", type=float, default=1., help="the gradient clipping bound")
    parser.add_argument("-s", "--sigma", type=float, default=1., help="the noise scale")
    parser.add_argument("-e0", "--eps0", type=float, default=4., help="the epsilon privacy budget")
    parser.add_argument("-d0", "--delta0", type=float, default=1e-05, help="the delta privacy budget")
    parser.add_argument("-S", "--save_paths", type=str, default=pm_paths, help="the path of the model to load")
    parser.add_argument("-G", "--G_path", type=str, default=pm_paths[0], help="the path of the generator model to load")
    parser.add_argument("-D", "--D_path", type=str, default=pm_paths[1], help="the path of the discriminator model to load")
    parser.add_argument("--tol", type=float, default=1e-08, help="the tolerance for training convergence")

    return parser.parse_args()


def pre_model(pretrained = False):
    im_shape = (78, 78, 1)
  
    p = pm_paths[0]
    if pretrained and os.path.exists(p):
        G = load_model(p)
    else:
        G = Sequential()
        units = np.logspace(*np.log2([10, 784]), 3, base=2).round().astype(int)
        G.add(Flatten(input_shape=im_shape))
        for u in reversed(units[:-1]):
            G.add(Dense(u, activation="relu", kernel_initializer="lecun_normal"))
        for u in units[1:]:
            G.add(Dense(u, activation="relu", kernel_initializer="lecun_normal"))
        G.add(Reshape(im_shape))
        G.compile("adam", "mse")

    p = pm_paths[1]
    if pretrained and os.path.exists(p):
        D = load_model(p)
    else:
        D = Sequential()
        D.add(Conv2D(64, 3, padding="same", activation="relu", kernel_initializer="lecun_normal", input_shape=im_shape))
        D.add(MaxPool2D(padding="same"))
        D.add(Conv2D(64, 3, padding="same", activation="relu", kernel_initializer="lecun_normal"))
        D.add(MaxPool2D(padding="same"))
        D.add(Conv2D(64, 3, padding="same", activation="relu", kernel_initializer="lecun_normal"))
        D.add(MaxPool2D(padding="same"))
        D.add(Flatten())
        D.add(Dense(512, activation="relu", kernel_initializer="lecun_normal"))
        D.add(Dense(512, activation="relu", kernel_initializer="lecun_normal"))
        D.add(Dense(10, activation="softmax"))
        D.compile("adam", "categorical_crossentropy", ["accuracy"])

    return G, D
    
    
    
def mean_estimation(y_true, y_pred):
    return K.mean(y_pred - y_true)

def evaluate(G, D, X, n = 10):
    D.compile('adam', mean_estimation)
    m = D.evaluate(X, np.ones((X.shape[0], 1)))
    print('Discriminator mean estimation of X_test: {:.4f}'.format(m))
    Z = np.random.uniform(size = X.shape)
    m = D.evaluate(Z, np.zeros((Z.shape[0],1)))
    print('Discriminator mean estimation of Z: {:.4f}'.format(m))
    Gz = G.predict(Z)
    m = D.evaluate(Gz, np.zeros((Gz.shape[0], 1)))
    print('Discriminator mean estimation of G(z): {:.4f}'.format(m))

    X_ = X[np.random.permutation(X.shape[0])[:n]]
    Z, Gz = Z[:n], Gz[:,n]
    y_pred = np.empty((3, n))
    y_pred[0, :] = D.predict_proba(X_).ravel()
    rho = np.random.uniform(size = n)
    X_hat = np.array([(rho[i]*X_[i]+(1-rho[i])*Z[i]).tolist() for i in range(n)])
    y_pred[1,:] = D.predict_proba(X_hat).ravel()
    y_pred[2,:] = D.predict_proba(Gz).ravel()

    fig, axs = plt.subplots(3, n, figsize = np.array((n, 3))*1., sharex = True, sharey = True)
    for i in range(3):
        for j in range(n):
            axis = axs[i, j]
            im = np.squeeze(X_[j] if i == 0 else X_hat[j] if i == 1 else Gz[j])
            axis.imshow(im, cmap = plt.cm.Greys)
            if j == 0:
                axis.set_title('D(x)' if i == 0 else 'D(\hat x)=' if i == 1 else 'D(G(z) =')
            axis.set_title("${}{:.2f}$".format(axis.get_title(), y_pred[i,j]))
        fig.tight_layout()
        return fig
            
            
            
if __name__ == "__main__":
    
    # load input
    results = []
    path = []
    path.append('D:/TUT/Medical/biophysics/experiment')
    with open('D:/TUT/Medical/biophysics/experiment/Data_Cortex_Nuclear.csv') as csvfile:
        reader = csv.reader(csvfile) # change contents to floats
        for row in reader: # each row is a list
            results.append(row)
    results = np.array(results)
    X = results[1:1081,0:78]
                                          

    args = parse_args()
    
    # split data
    X_train, X_test = train_test_split(X,test_size=0.2)
    
    if (os.path.exists(pm_paths[0]) & os.path.exists(pm_paths[1])):
        G, D = pre_model(pretrained = True)
    else:
        G, D = pre_model(pretrained = False)

    if os.path.exists(args.G_path):
        print('Loading model G')
        G = load_model(args.G_path)
    if os.path.exists(args.D_path):
        print('Loading model D')
        D = load_model(args.D_path)
    
    # Adapting the discriminator net
    D.add(Dense(1, activation = 'sigmoid', name = 'real'))    
    # for i in range(len(D.layers)):
    #    if isinstance(D.layers[i], Conv2D):
    #           D.layers[i].trainable = False
    
    try:
        # dp-GAN training
        if args.epochs > 0:
            save_paths = []
            for path, net in zip([args.G_path, args.D_path], ['G', 'D']):
                dir, base = os.path.split(path)
                if 'DpGAN' not in base:
                    base = '_'.join(('DpGAN', net))
                save_paths.append(os.path.join(dir, base))
            train_kw = {'epochs': args.epochs,
                        'sess': args.sess,
                        'zdim': args.zdim,
                        'img_size': args.img_size,
                        'lam': args.lam,
                        "n_critic": args.critic,
                        "batch_size": args.batch,
                        "optimizer": tf.train.AdamOptimizer(*args.adam),
                        "C": args.clip,
                        "sigma": args.sigma,
                        "eps0": args.eps0,
                        "delta0": args.delta0,
                        "save_paths": args.save_paths,
                        "tol": args.tol}
 
#            dpgan = DpGAN(G, D)
#            G, D, G_loss_hist, D_loss_hist = dpgan.train(X_train, **train_kw, out_path)     
            G, D, G_loss_hist, D_loss_hist, t_train = run(X_train)

    finally:
        print("Evaluating the nets")
        fig = evaluate(G, D, X_test)
        path = "summary/dpgan_evaluation_t{}_l{}_a{}_b1{}_b2{}_C{}_s{}".format(args.critic, args.lam, *args.adam, args.clip, args.sigma)
        if os.path.exists(path + ".png"):
            path += "_bis"
        path += ".png"
        fig.savefig(path)

        # Plotting the training losses
        fig, ax = plt.subplots()
        xx = np.arange(len(D_loss_hist))+1
        ax.plot(xx[::args.critic]+args.critic-1, G_loss_hist, "-o", label="G's loss", ms=3.6)
        ax.plot(xx, D_loss_hist, "-x", label="D's loss", ms=3.6)
        ax.grid(axis="y")
        ax.set_xlabel("Batch iteration")
        ax.set_title(r"Parameters: $n_\mathit{critic}"+r"={}, \lambda={}, C={}, \sigma={},$".format(args.critic, args.lam, args.clip, args.sigma)+"\n"+r"$\alpha={}, \beta_1={}, \beta_2={}$".format(*args.adam))
        ax.legend(loc="upper right")
        path = "summary/dpgan_loss_history_t{}_l{}_a{}_b1{}_b2{}_C{}_s{}".format(args.critic, args.lam, *args.adam, args.clip, args.sigma)
        if os.path.exists(path + ".png"):
            path += "_bis"
        path += ".png"
        fig.savefig(path)

        print("Close the figures to quit")
        plt.show()
            
