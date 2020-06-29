# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 09:38:21 2019

@author: user
"""
import abc
from datetime import datetime as dt
import numpy as np
import tensorflow as tf
from six.moves import xrange
from test import *
#from keras.datasets import mnist


# This section is to implement DCGAN. You need to complete TODO parts by yourself according to DCGAN topology
class BasicDpGAN(abc.ABC):
    def _init_(self, G, D):
        self.G = G
        self.D = D


    @abc.abstractclassmethod
    def train(self):
        pass

    def save(self, paths):
        self.G.save(paths[0])
        self.D.save(paths[1])

class DpGAN(BasicDpGAN):
    def __init__(self, G, D, sess, img_size, n_critic, z_dim, lam, batch_size, epochs, C, sigma, eps0, delta0, tol, out_paths):
        super().__init__(G, D)
            
        self.sess = sess
        self.epochs = epochs
        self.z_dim = z_dim
        self.img_size = img_size
        self.img_dim = img_size * img_size
        self.img_shape = [img_size, img_size, 1]
        self.batch_size = batch_size
        self.n_critic = n_critic,
        self.lam = lam
        self.C = C
        self.sigma = sigma, 
        self.eps0 = eps, 
        self.delta0 = delta0
        self.tol = tol
        self.out_paths = out_paths
        self.build_model()
        self.model_name = "DpGAN.model"

    def build_model(self):
        self.is_training = tf.placeholder(tf.bool, name="is_training")
        self.img = tf.placeholder(tf.float32, [None]+self.img_shape, name='real_images')

        self.z = tf.placeholder(tf.float32, [None, self.z_dim], name='z')
        self.z_sum = tf.summary.histogram('z', self.z)

#        initializer = tf.contrib.layers.xavier_initializer()

        # self.img_fake is produced by generator with a random input z.
        self.rho, self.img_fake = self.generator(self.z)
        self.x_hat = tf.multiply(self.rho, self.img) + tf.multiply(tf.constant(1.) - self.rho, self.G(self.img_fake))
        # Outputs from discriminator with real image or fake image.
        self.D_real, self.D_logits_real = self.discriminator(self.img)
        self.D_fake, self.D_logits_fake = self.discriminator(self.img_fake)

        self.d_sum = tf.summary.histogram("d", self.D)
        self.d__sum = tf.summary.histogram("d_", self.D_fake)
        self.img_fake_sum = tf.summary.histogram("G", self.img_fake)

        t_vars = tf.trainable_variables()
        self.d_vars = [var for var in t_vars if var.name.startswith('discriminator')]
        self.g_vars = [var for var in t_vars if var.name.startswith('generator')]

        # Calculating Loss value
        self.D_grad_x_hat = tf.gradients(self.D(self.x_hat))[0] 
#        self.D_loss_real = tf.reduce_mean(self.D(self.G(self.D_logits_real) - self.rho + tf.multiply(tf.constant(tf.float32, lam), tf.square(tf.norm(D_grad_x_hat, axis = 1 if x_hat._rank() < = 2 else(1, 2)) - tf.constant(1.))), axis = 0)
#        self.D_loss_fake = tf.reduce_mean(self.D(self.G(self.D_logits_fake) - self.rho + tf.multiply(tf.constant(tf.float32, lam), tf.square(tf.norm(D_grad_x_hat, axis = 1 if x_hat._rank() < = 2 else(1, 2)) - tf.constant(1.))), axis = 0)
        self.D_loss_real = tf.reduce_mean(self.D(self.G(self.D_logits_real) - self.rho)) 
        self.D_loss_fake = tf.reduce_mean(self.D(self.G(self.D_logits_fake) - self.rho)) 
        self.D_loss = self.D_loss_real + self.D_loss_fake
#        self.D_grad_vars = optimizer.compute_gradients(self.D_loss, self.D.trainable_weights)
        #clipping and noising
#        self.D_grad_vars = [(tf.divide(grad, tf.maximum(tf.constant(1.), tf.divide(tf.divide(tf.norm(grad), tf.constant(tf.float32, C)))))+ tf.random_normal(grad.shape, 0, sigma * C), var) for grad, var in D_grad_vars] 
#        self.train_D = optimizer.apply_gradients(self.D_grad_vars)

        self.G_loss = tf.reduce_mean(-self.D(self.G(self.z)), axis = 0)
#        self.train_G = self.optimizer.minimize(self.G_loss, var_list = self.G.trainable_weights)
        
    def generator(self, z, batch_size,):

        #generate the random input as size z
        rho = tf.random_uniform(tf.float32, (self.bach_size,)) 
        for d in z.shape[1,:].as_list(): rho = tf.stack(d * [z], axis = -1)
        logits = tf.random_uniform((batch_size,) + self.G.input_shape[1:], name = 'z')
        return rho, logits

    def discriminator(self, img, batch_size,):

        #generate the random input as size z
        rho = tf.random_uniform(tf.float32, (self.bach_size,)) 
        for d in img.shape[1,:].as_list(): rho = tf.stack(d * [img], axis = -1)
        logits = tf.random_uniform((batch_size,) + self.G.input_shape[1:], name = 'z')
        return rho, logits
    

    def train(self, X, save_paths):
        # initialize the variables
        delta = 0
        d_optim = tf.train.AdamOptimizer().minimize(self.D_loss, X, var_list=self.d_vars)
        g_optim = tf.train.AdamOptimizer().minimize(self.G_loss, X, var_list=self.g_vars)

        try:
            tf.global_variables_initializer().run()
            tf.local_variables_initializer().run()
        except:
            tf.initialize_all_variables().run()
#        batch_loss_hist = []
        D_loss_hist = []
        G_loss_hist = []
        S = np.shape(X)
        N = S[0]
        for n in xrange(N+1):
            Xn = X[n, :]            
            for epoch in xrange(self.epochs+1):
                input_z = np.random.uniform(-1, 1, [self.batch_size, self.z_dim]).astype(np.float32)
                input_imgs, _ = Xn.train.next_batch(self.batch_size)
                input_imgs_ = input_imgs.reshape(self.batch_size, self.img_size, self.img_size, 1)
                for t in xrange():
                # Update D network and G network
                    weights = tf.get_variable("weights", shape=(), initializer=tf.zeros_initializer())
                    weights = weights.assign_add(weights)
                    weights.assign(weights)
                    theta_old = self.sess.run([weights for weights in self.G.trainable_weights])
                    _, [D_loss_curr] = self.sess.run([d_optim, self.D_loss], feed_dict={self.img: input_imgs_, self.z: input_z})
                    _, [G_loss_curr] = self.sess.run([g_optim, self.G_loss], feed_dict={self.z: input_z})
                    D_loss_hist.append(D_loss_curr)            
                    G_loss_hist.append(G_loss_curr)            
                    
                t_vars = tf.trainable_variables()
                d_vars = [var for var in t_vars if var.name.startswith('discriminator')]
                g_vars = [var for var in t_vars if var.name.startswith('generator')]
    
                #save model
                if save_paths:
                    self.save(save_paths)
                
                # Stopping criterias
                if self.tol >= 0:
                    G_convergence = all(self.sess.run(tf.norm(self.G.trainable_weights[i] - theta_old[i]) <= self.tol) for i in range(len(self.G.trainable_weights)))
                    if (delta > self.delta0) or G_convergence:
                        print("Epoch: [{:2d}/{:2d}] -D_loss: {:.8f}, -G_loss {:.8f} /n".format(epoch, self.epochs, D_loss_curr, G_loss_curr))
                        break
    #              samples = self.sess.run(self.img_fake, feed_dict={self.z: input_z})
            # Returning results
            if self.tol >= 0:
                if G_convergence:
                    print('Training stopped: G has converged')
                if delta > self.delta0:
                    print('Training stopped: Privacy budget exceeded (delta = {})'.format(delta))
        self.sess.close()
        return self.G, self.D, G_loss_hist, D_loss_hist
                #import pdb; pdb.set_trace()
#                fig = self.plot(samples)
#                if not os.path.exists(output_path):
#                    os.makedirs(output_path)
#                plt.savefig(os.path.join(output_path, '{}.png'.format(str(i).zfill(3))), bbox_inches='tight')
#                i += 1
#                plt.close(fig)

         
def run(X):
    config = tf.ConfigProto(device_count = {'GPU': 0}) # If you wanna implement this code with GPU, change 0 to 1 or 2
    # pre-training
#    dataset = input_data.read_data_sets('./MNIST_data', one_hot=True)

    with tf.Session(config=config) as sess:
        # dcGAN
        tic = dt.now()
        datasize = np.shape(X)
        try:
            dpgan = DpGAN(args.G, args.D, sess = args.sess, img_size=datasize[1], n_critic = args.n_critic, z_dim=datasize[0], lam = args.lam, batch_size=args.batch_size, epochs=args.epochs, C = args.C, sigma = args.sigma, eps0 = args.eps0, delta0 = args.delta0, tol = args.tol, out_paths = args.save_paths)
            G, D, G_loss_hist, D_loss_hist = dpgan.train(X, args.save_path)
        except Exception:
            print('arguments error')
        toc = dt.now() - tic

    return G, D, G_loss_hist, D_loss_hist, toc


