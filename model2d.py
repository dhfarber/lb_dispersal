# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 12:10:53 2018

@author: daniel.farber
"""
import numpy as np
from math import sqrt,exp
from scipy.special import gamma
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d as mp3

class kernel:
    def __init__(self,n):
        self.n = n
    def initialize(self):
        return(np.ones(self.n))
    def inv_power(self,b):
        k = np.zeros(self.n)
        for i in range(self.n):
            k[i] = (i+0.5)**(-b) 
        y = k/(k.sum())
        return(y)
    def gamma_dist(self,a,b): 
        '''issues at dist = 0'''
        k = np.zeros(self.n)
        for i in range(self.n):
            k[i] = ((b**a)/gamma(a))*(i**(a-1))*exp(-b*i)
        y = k/(k.sum())
        return(y[1:])  
def dist_from_source(source_x,source_y,field_x,field_y): #returns array of approximate distance each cell in an arraty 
    cool = np.ones((2*field_x,2*field_y))
    return_arr = np.zeros((2*field_x,2*field_y))
    for i, j in np.argwhere(cool):
        return_arr[i,j] = round(sqrt((source_x + round((field_x/2)) - i)**2 + (source_y + round((field_y/2))- j)**2),0)
    return(return_arr.astype(int))
def kernel_2d(dist_arr,x,y,kernel_1d): 
    return_arr = np.ones_like(dist_arr); #get double the size of the field to account for spores blown out of the field
    pdf = kernel_1d
    for i,j in np.argwhere(return_arr):
        return_arr[i,j] = pdf[int(dist_arr[i,j])] #return_arr[i,j] = pdf[int(dists[int(round(i+(x/2))),int(round(j+(y/2)))])] #currently returns correct size and shape array, but doesn't take into account spores traveling outside of field
    normed_ret_arr = return_arr/np.sum(return_arr) #normalize each cell so that entire 2x by 2y array adds up to 1
    return(normed_ret_arr[round(.1+(x)/2):round(.1+(x*(3/2))),round(.1+(y/2)):round(.1+(y*(3/2)))])
def model_2d(x,y,T,spacing,l_period,i_period,dmfr,sites,initial_sev,focus,distribution,parameters): 
    if distribution=='inverse power' or distribution=='inv_power':
        kernel_1d = kernel(3*x).inv_power(parameters) #one-dimensional kernel
    elif distribution=='gamma':
        kernel_1d = kernel(3*x).gamma_dist(parameters[0],parameters[1]) #one-dimensional kernel
    else:
        return('Error: kernel not recognized. Please choose "gamma" or "inverse power"')
    dists_arr = np.empty((2*x,2*y,x*y)) #array of distances from sources in x by y field to 2x by 2y potential "sinks" 
    counter = 0
    for i,j in np.argwhere(np.ones((x,y))): 
        dists_arr[:,:,counter] = dist_from_source(i,j,x,y) 
        counter += 1
    kern2d_arr = np.empty((x,y,x*y)) #2D dispersal kernel for all possibile source/sink combinations. 
    kern2d_arr[:,:,0] = kernel_2d(dists_arr[:,:,0],x,y,kernel_1d)
    for i in range(0,x*y):
        kern2d_arr[:,:,i] = kernel_2d(dists_arr[:,:,i],x,y,kernel_1d)
    p_return = np.zeros((T+1,x,y,1)) #healthy array p empty array
    for i in range(0,x,spacing): p_return[:,:,i] = sites 
    p_return[0,focus[0],focus[1]] = sites - sites*initial_sev #p with healthy plants - initial latent infection
    u_return =  np.zeros((T+1,x,y,l_period)) #latent array u empty
    u_return[0,focus[0],focus[1],0] = sites*initial_sev #latent array day 0
    v_return =  np.zeros((T+1,x,y,i_period)) #infectious array v - all 0s (this is the state of v on day 0)
    rem_return =  np.zeros((T+1,x,y,1)) #removed array rem - all 0s (this is the state of rem on day 0)
    for t in range(T):
        print('day ' + str(t+1))
        ret_arr = u_return[t,:,:,:] #copy latent array from previous day T to start
        ret_arr[:,:,1:] = u_return[t,:,:,0:l_period-1]
        lat_sum = np.sum(u_return[t,:,:,:],axis=2)
        inf_sum = np.sum(v_return[t,:,:,:],axis = 2)
        infec_sum = inf_sum #copy of inf_sum to manipulate for creating new latent lesions
        dis = sum(lat_sum,inf_sum)
        p_return[t,:,:] = p_return[t,:,:] - dis.reshape(x,y,1) - rem_return[t,:,:] #healthy faction for t
        p_return[p_return<0] = 0
        for i,j in np.argwhere(np.ones_like(ret_arr[:,:,0])):
            new_inf = (dmfr*infec_sum*kern2d_arr[:,:,int(x*i+j)])/sites #this is it.
            sum_new_inf = np.sum(new_inf)
            ret_arr[i,j,0] = min(1,sum_new_inf)*p_return[t,i,j]
        u_return[t+1,:,:,:] = ret_arr
        i_new_les = np.zeros_like(v_return[t,:,:,:])
        i_new_les[:,:,0] = u_return[t,:,:,l_period-1].reshape(x,y)
        i_new_les[:,:,1:] = v_return[t,:,:,0:i_period-1]
        v_return[t+1:,:,:] = i_new_les
        rem_return[t+1,:,:] = v_return[t,:,:,i_period-1].reshape(x,y,1)
        print('Mean healthy tissue remaining:')
        print(round(np.mean(p_return[t,:,:]),2))
        if(np.mean(p_return[t,:,:])<.03*sites): #break when disease severity is > 97%
            break
    rem_return = rem_return.reshape(T+1,x,y)
    ret_dir = {'p':(p_return),'u':(u_return),'v':(v_return),'rem':( rem_return),'k':(kernel_1d)}
    return(ret_dir)
class plot_model2d:
    def __init__(self,ret_dir): #input return dictionary from model_2d_start
        self.ret_dir = ret_dir
        self.dis = np.sum(ret_dir['u'],axis=3)+np.sum(ret_dir['v'],axis=3)+ret_dir['rem'].reshape(len(ret_dir['rem'][:,1,1]),len(ret_dir['rem'][1,:,1]),len(ret_dir['rem'][1,1,:]))
        #print(self.dis)
    def plot_1d(self,latent,focus): #latent is a scalar; neceessary to plot a line at the end of each latent period
        #calculate disease as a functgion of distance?
        dist_arr = dist_from_source(focus[0],focus[1],np.size(self.dis,axis = 0),np.size(self.dis, axis = 1))
        plotting_arr = np.zeros((np.amax(dist_arr),np.size(self.dis,axis=2)))
        #print(plotting_arr)
        for t in range(np.size(self.dis,axis=2)):#range(latent,np.size(self.dis,axis=2), latent): #issue: going from 2d to 1d
            for d in range(np.amax(dist_arr)):
                #print(d)
                for i, j in np.argwhere(dist_arr==d): #pull row and column indices for each distance 
                   plotting_arr[:,t] = sum(self.dis[i,j])
        nmb_latent_per = round(len(plotting_arr[1,:])/latent)
        col_list = ['purple','navy','blue','c','green','yellow','orange','red','brown','black']
        for i in range(nmb_latent_per):
            plt.plot(range(np.amax(dist_arr)),plotting_arr[:,latent*i],color = col_list[i])
        return(plotting_arr)
    def plot_2d(self,x,y,day):#plots 3d figure of x = N-S, y = W-E, z = severity at given day
        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')
        X,Y = np.meshgrid(np.linspace(0,x,num=x),np.linspace(0,x,num=x))        
        ax.plot_surface(X,Y,self.dis[day,:,:],cmap = "coolwarm")
        title = "Day: " + str(day)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes,fontsize = 18)
        #ax.plot_surface(range(x),range(y),self.dis[day,:,:],cmap = "coolwarm")
        return(self.dis[day,:,:])
