# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 12:10:53 2018

@author: daniel.farber
"""
'''code optimized by get kernel prior to iterating through days of epidemic. 09-26-18'''
'''currently fields must be square. 09-28-18'''
'''lesion growth works as expected. 09-28-18'''
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
    def inv_power(self,b,c):
        k = np.zeros(self.n)
        for i in range(self.n):
            k[i] = (i+c)**(-b) 
        y = k/(k.sum())
        return(y)
    def gamma_dist(self,a,b): 
        '''issues at dist = 0'''
        k = np.zeros(self.n)
        for i in range(self.n):
            k[i] = ((b**a)/gamma(a))*((i+.0000000001)**(a-1))*exp(-b*(i+.0000000001)) #trial fix: adding tiny amount to each distance observation
        y = k/(k.sum())
        return(y)  
def dist_from_source(source_x,source_y,field_x,field_y): #returns array of approximate distance of each cell in an array 
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
def model_2d(x,y,T,spacing,l_period,i_period,dmfr,sites,initial_sev,focus,distribution,parameters,lesion_growth): 
    if distribution=='inverse power' or distribution=='inv_power':
        kernel_1d = kernel(3*x).inv_power(parameters[0],parameters[1]) #one-dimensional kernel
    elif distribution=='gamma':
        kernel_1d = kernel(3*x).gamma_dist(parameters[0],parameters[1]) #one-dimensional kernel
    else:
        return('Error: kernel not recognized. Please choose "gamma" or "inverse power"')
    dists_arr = np.empty((2*x,2*y,x*y)) #array of distances from sources in x by y field to 2x by 2y potential "sinks"
    dists_field_subset = dists_arr[round(.1+(x)/2):round(.1+(x*(3/2))),round(.1+(y/2)):round(.1+(y*(3/2)))]
    counter = 0
    for i,j in np.argwhere(np.ones((x,y))): 
        dists_arr[:,:,counter] = dist_from_source(i,j,x,y) 
        counter += 1
    kern2d_arr = np.empty((x,y,x*y)) #2D dispersal kernel for all possibile source/sink combinations. 
    kern2d_arr[:,:,0] = kernel_2d(dists_arr[:,:,0],x,y,kernel_1d)
    for i in range(0,x*y):
        kern2d_arr[:,:,i] = kernel_2d(dists_arr[:,:,i],x,y,kernel_1d)
    p_return = np.zeros((T+1,x,y,1)) #healthy array p empty array
    for i in range(0,x,spacing): p_return[:,i,:] = sites 
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
        #print('latent: ' + str(round(np.mean(lat_sum),2)))
        inf_sum = np.sum(v_return[t,:,:,:],axis = 2)
        #print('infectious: ' + str(round(np.mean(inf_sum),2)))
        infec_sum = inf_sum #copy of inf_sum to manipulate for creating new latent lesions
        dis = sum(lat_sum,inf_sum)
        p_return[t,:,:] = p_return[t,:,:] - dis.reshape(x,y,1) - rem_return[t,:,:] #healthy faction for t
        p_return[p_return<0] = 0
        #print('Mean healthy tissue remaining:')
        #print(round(np.mean(p_return[t,:,:]),2)) #there's a small error here somewhere. healthy tissue goes up on final day of latent period?
        i_new_les = np.zeros_like(v_return[t,:,:,:]) #infectious lesions empty array
        i_new_les[:,:,0] = u_return[t,:,:,l_period-1].reshape(x,y) #new infectious lesions
        i_new_les[:,:,1:] = v_return[t,:,:,0:i_period-1] #infectious lesion growth
        for i,j in np.argwhere(np.ones_like(ret_arr[:,:,0])):
            new_inf = (dmfr*infec_sum*kern2d_arr[:,:,int(x*j+i)])/sites #this is it.
            sum_new_inf = np.sum(new_inf) #sum of total infections (as a fraction of 1)
            ret_arr[i,j,0] = min(1,sum_new_inf)*p_return[t,i,j] #multiplies the healthy tisse remaining by the the fraction to be newly infected
            for d in range(1,i_period): #for each day of infectious period, add lesion_growth to each lesion, provided that it doesn't add up to more than total leaf area
                if np.sum(i_new_les[i,j,:]) + lesion_growth <  p_return[t,i,j] - np.sum(ret_arr[i,j]) - rem_return[t,i,j] and i_new_les[i,j,d] >= 1: #cutoff at 1. will work better stochastically
                    i_new_les[i,j,d] += lesion_growth
        u_return[t+1,:,:,:] = ret_arr
        v_return[t+1:,:,:] = i_new_les
        rem_return[t+1,:,:] = rem_return[t,:,:] + v_return[t,:,:,i_period-1].reshape(x,y,1)    #fixed rem_return
        #if(np.mean(p_return[t,:,:])<.02*sites): #break when disease severity is > 97%
         #   break
    rem_return = rem_return.reshape(T+1,x,y)
    p_return = p_return.reshape(T+1,x,y)
    l_sum = np.sum(u_return,axis = 3)
    v_sum = np.sum(v_return,axis = 3)
    diseased = l_sum + v_sum + rem_return
    ret_dir = {'Diseased': diseased, 'Healthy':(p_return),'Latent':(u_return),'Latent sum':l_sum,'Infectious':(v_return),'Infectious Sum':v_sum,'Removed':( rem_return),'kernel':(kernel_1d),'distances':(dists_field_subset),'x':x,'y':y,'T':T,'spacing':spacing,'latent period':l_period,"infectious period":i_period,'DMFR':dmfr,'sites':sites,'initial severity':initial_sev,'focus':focus,'distribution':distribution,'paramters':parameters,'lesion growth':lesion_growth}
    return(ret_dir)
class plot_model2d:
    def __init__(self,ret_dir): #input return dictionary from model_2d_start
        self.ret_dir = ret_dir
        self.dis = np.sum(ret_dir['Latent'],axis=3)+np.sum(ret_dir['Infectious'],axis=3)+ret_dir['Removed'].reshape(len(ret_dir['Removed'][:,1,1]),len(ret_dir['Removed'][1,:,1]),len(ret_dir['Removed'][1,1,:]))
        self.distance_from_focus = self.ret_dir['distances'][:,:,int(self.ret_dir['focus'][0] + (self.ret_dir['focus'][1] * self.ret_dir['x']))].astype(int) # get 2D array of distansces from initial focus to all other cells in field
    def plot_1d(self): #plot the downwind cells self.dis as a function of the number of cells
        col_list = ['purple','navy','blue','c','green','yellow','orange','red','brown','black'] #colors to represent the disease progress at each latent period
        focus_row = self.ret_dir['focus'][1] 
        x_arr = np.arange(focus_row,len(self.dis[0,:,:])) # 1D x array
        nmb_latent_per = round(float(self.ret_dir['T']/self.ret_dir['latent period'])) #number of latent periods
        for i in range(nmb_latent_per):
            disease_downwind = self.dis[:,focus_row, self.ret_dir['focus'][0]:] #1D y array
            plt.plot(x_arr, disease_downwind[self.ret_dir['latent period']*i,:],color = col_list[i])
        plt.title('Disease gradient from initial focus')
        plt.xlabel('Distance from source')
        plt.ylabel('Number of infections')
        plt.legend('Latent periods' + range(nmb_latent_per))
        #return(disease_downwind)
    def plot_2d(self,x,y,day):#plots 3d figure of x = N-S, y = W-E, z = severity at given day
        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')
        X,Y = np.meshgrid(np.linspace(0,x,num=x),np.linspace(0,x,num=x))        
        ax.plot_surface(X,Y,self.dis[day,:,:],cmap = "coolwarm")
        title = "Day: " + str(day)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes,fontsize = 18)
        #ax.plot_surface(range(x),range(y),self.dis[day,:,:],cmap = "coolwarm")
        return(self.dis[day,:,:])