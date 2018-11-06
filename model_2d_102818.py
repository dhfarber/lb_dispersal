# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 18:29:12 2018

@author: Daniel
"""
'''102718 - uses kern_2d which gets max(x,y)**2 array of which (x,y) has contains > 0 sites'''
'''code optimized by get kernel prior to iterating through days of epidemic. 09-26-18'''
'''currently fields must be square. 09-28-18'''
'''lesion growth works as expected. 09-28-18'''
'''something is wrong when p is very close to 0, so that both the infectious array and the latent array produce a sinosoidal pattern'''
'''trouble-shooting the saturation issue: the day when p should reach 0, it instead pops way up'''
'''also, weird thing: whole row where initial focus is is reduced'''
'''thought: p[t+1] = sites - (u[t+1] + v[t+1] + rem[t+1]) because it should be based on the present day, not based on some function of p[t]'''
import numpy as np
from math import sqrt,exp
from scipy.special import gamma
import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d as mp3

class kernel:
    '''creates a 1D PDF object'''
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

def dist_from_source(source_x, source_y, field):  
    '''returns array of approximate distance of each cell in a
    [2 * max(x, y), 2 * max(x, y)] array from a source located 
    in cell [source_x, source_y]'''
    y = len(field[:, 0])
    x = len(field[0, :])
    max_length = max(x, y)
    dummy = np.ones((max_length * 2, max_length * 2))  # get square array with sides = max length of field
    dists_arr = np.zeros_like(dummy)
    # return_arr should return distances from source(x,y) in field to each cell in return_arr
    # this means centering field
    # therefore both source_x and source_y are displaced by half of the length of the longest dimension of field
    for i, j in np.argwhere(dummy):
        dists_arr[i, j] = round(sqrt((source_x + 
                 round(((max_length + .1) / 2)) - i) ** 2 + 
                 (source_y + round(((max_length + .1) / 2)) - j) ** 2), 0)
    return dists_arr.astype(int)

def kernel_2d(dist_arr, x, y, kernel_1d):
    '''returns an array the size (x,y) of fraction of total new lesions in 
    each cell. dist_arr is a square array of length 2*max(x,y) on each side
    and y are inputs from model for the simulated field kernel_1d is the 
    one-dimensional dispersal kernel as a function of distance from a source'''
    field_arr = np.ones_like(dist_arr)
    return_arr = np.ones_like(field_arr,float)
    pdf = kernel_1d
    for i, j in np.argwhere(field_arr):
        return_arr[i, j] = pdf[int(dist_arr[i, j])]
    normed_ret_arr = return_arr/np.sum(return_arr) #normalize each cell so that entire 2x by 2y array adds up to 1
    max_length = max(x,y)
    kern2d_sq = np.zeros((max_length,max_length))
    kern2d = normed_ret_arr[round(.1 + max_length / 2):round(.1 + max_length / 2) + y,
                        round(.1 + max_length / 2):round(.1 + max_length / 2) + x]
    kern2d_sq[:y,:x] = kern2d
    return kern2d_sq

def kernel_array(kernel_1d, field):
    '''returns 2 3D arrays: 1) of all 2D distance_from_srouce's
    2) all 2D kernels from all possible sources.'''
    y = len(field[:, 0])
    x = len(field[0, :])
    max_length = max(x, y)
    dists_arr = np.ones((2 * max_length,2 * max_length, max_length**2))
    counter = 0
    field = np.ones((max_length, max_length))
    for i,j in np.argwhere(field):
        dists_arr[:,:,counter] = dist_from_source(i,j,field)
        counter += 1
    kern2d_arr = np.zeros((max_length, max_length, max_length**2))
    for i in range(0, x * y ):
        kern2d_arr[:, :, i] =  kernel_2d(dists_arr[:, :, i], 
                  x , y, kernel_1d)
    return dists_arr, kern2d_arr
def model_2d(x,y,T,spacing,l_period,i_period,dmfr,sites,initial_sev,focus,distribution,parameters,lesion_growth): 
    '''models spread of an invading population'''
    max_length = max(x, y)
    if distribution=='inverse power' or distribution=='inv_power':
        kernel_1d = kernel(3*x).inv_power(parameters[0],parameters[1]) #one-dimensional kernel
    elif distribution=='gamma':
        kernel_1d = kernel(3*x).gamma_dist(parameters[0],parameters[1]) #one-dimensional kernel
    else:
        return('Error: kernel not recognized. Please choose "gamma" or "inverse power"')
    dists3d, kern3d = kernel_array(kernel_1d, np.ones((x, y)))
    dists_field_subset = dists3d[round(.1 + max_length / 2):round(.1 + max_length / 2) + y,
                                 round(.1 + max_length / 2):round(.1 + max_length / 2 + x)]
    p_return = np.zeros((T+1,max_length,max_length,1)) #healthy array p empty array
    for i in range(0,x,spacing): p_return[:y,i,:] = sites
    #p_return[0,focus[0],focus[1]] = sites - sites*initial_sev #p with healthy plants - initial latent infection
    u_return =  np.zeros((T+1,max_length,max_length,l_period+1)) #latent array u empty
    u_return[0,focus[0],focus[1],0] = sites*initial_sev #latent array day 0
    v_return =  np.zeros((T+1,max_length,max_length,i_period+1)) #infectious array v - all 0s (this is the state of v on day 0)
    rem_return =  np.zeros((T+1,max_length,max_length,1)) #removed array rem - all 0s (this is the state of rem on day 0)
    for t in range(T):
        print('day ' + str(t+1))
        dis = np.sum(u_return[t,:,:,:], axis = 2) + np.sum(v_return[t,:,:,:], axis = 2)
        p_return[t,:,:] = sites - dis.reshape(max_length,max_length,1) - rem_return[t,:,:] #healthy faction for t
        p_return[p_return<0] = 0
        for l in range(1, l_period+1):
            u_return[t+1,:,:,l] = u_return[t,:,:,l-1]
        inf_sum = np.sum(v_return[t,:,:,:],axis = 2)
        infec_sum = inf_sum #copy of inf_sum to manipulate for creating new latent lesions
        i_new_les = np.zeros_like(v_return[t,:,:,:]) #infectious lesions empty array
        i_new_les[:,:,1:] = v_return[t,:,:,0:i_period] #infectious lesion growth
        i_new_les[:,:,0] = u_return[t,:,:,l_period].reshape(max_length,max_length) #new infectious lesions
        for i,j in np.argwhere(np.ones_like(u_return[0,:,:,0])):
            new_inf = (dmfr * infec_sum * kern3d[:,:,int(x*i+j)])/sites 
            sum_new_inf = np.sum(new_inf) #sum of total infections (as a fraction of 1)
            if np.any(new_inf > 1):
                print(new_inf[new_inf>1])
            u_return[t+1,i,j,0] = min(1,sum_new_inf) * p_return[t,i,j] #multiplies the healthy tisse remaining by the the fraction to be newly infected
            for d in range(1,i_period): 
                if np.sum(i_new_les[i,j,:]) + lesion_growth <  p_return[t,i,j] - np.sum(u_return[t,i,j,:]) - rem_return[t,i,j] and i_new_les[i,j,d] >= 1: #cutoff at 1. will work better stochastically
                    i_new_les[i,j,d] += lesion_growth
        v_return[t+1:,:,:] = i_new_les
        rem_return[t+1,:,:] = rem_return[t,:,:] + v_return[t,:,:,i_period-1].reshape(max_length, max_length,1)    #fixed rem_return
        p_return[t,:,:] = sites - dis.reshape(max_length,max_length,1) - rem_return[t,:,:] #healthy faction for t
        p_return[p_return<0] = 0
    rem_return = rem_return.reshape(T+1, max_length, max_length)
    p_return = p_return.reshape(T+1, max_length, max_length)
    l_sum = np.sum(u_return,axis = 3)
    v_sum = np.sum(v_return,axis = 3)
    diseased = l_sum + v_sum + rem_return
    ret_dir = {'Diseased': diseased, 'Healthy':p_return,'Latent':u_return,'Latent sum':l_sum,'Infectious':v_return,'Infectious Sum':v_sum,'Removed':rem_return,'Kernel':kernel_1d,'Kernel Array':kern3d, 'Distances':dists_field_subset,'x':x,'y':y,'T':T,'Spacing':spacing,'Latent Period':l_period,"Infectious Period":i_period,'DMFR':dmfr,'Sites':sites,'Initial Severity':initial_sev,'Focus':focus,'Distribution':distribution,'Paramters':parameters,'Lesion Growth':lesion_growth}
    return(ret_dir)
class plot_model2d:
    def __init__(self,ret_dir): #input return dictionary from model_2d_start
        self.ret_dir = ret_dir
        self.dis = np.sum(ret_dir['Latent'],axis=3)+np.sum(ret_dir['Infectious'],axis=3)+ret_dir['Removed'].reshape(len(ret_dir['Removed'][:,1,1]),len(ret_dir['Removed'][1,:,1]),len(ret_dir['Removed'][1,1,:]))
        self.distance_from_focus = self.ret_dir['Distances'][:,:,int(self.ret_dir['Focus'][0] + (self.ret_dir['Focus'][1] * self.ret_dir['x']))].astype(int) # get 2D array of distansces from initial focus to all other cells in field
    def plot_1d(self,l_pers_plotted): #plot the downwind cells self.dis as a function of the number of cells
        col_list = ['purple','navy','blue','c','green','yellow','orange','red','brown','black'] #colors to represent the disease progress at each latent period
        focus_row = self.ret_dir['Focus'][0] 
        print('focus row: ' + str(focus_row))
        x_arr = np.arange(focus_row, len(self.dis[0,:,:])) # 1D x array
        print(x_arr)
        #nmb_latent_per = round(float(self.ret_dir['T']/self.ret_dir['Latent Period'])) #number of latent periods
        #for i in range(nmb_latent_per):
        for i in range(l_pers_plotted):
            disease_downwind = self.dis[:,focus_row, self.ret_dir['Focus'][0]:] #1D y array
            plt.plot(x_arr, disease_downwind[self.ret_dir['Latent Period']*i,:],color = col_list[i])
        plt.title('Disease gradient. DMFR: ' + str(self.ret_dir['DMFR']) + 
                  '; Initial severity: ' + str(self.ret_dir['Initial Severity']) + 
                  '%; b: ' + str(self.ret_dir['Paramters'][0]) + '; c: ' + 
                  str(self.ret_dir['Paramters'][1])) 
        plt.xlabel('Distance from source')
        plt.ylabel('Number of infections')
        return(disease_downwind)
    def plot_2d(self,x,y,day):#plots 3d figure of x = N-S, y = W-E, z = severity at given day
        fig = plt.figure(figsize=(10,10))
        ax = fig.gca(projection='3d')
        X,Y = np.meshgrid(np.linspace(0,x,num=x),np.linspace(0,x,num=x))        
        ax.plot_surface(X,Y,self.dis[day,:,:],cmap = "coolwarm")
        title = "Day: " + str(day)
        ax.text2D(0.05, 0.95, title, transform=ax.transAxes,fontsize = 18)
        #ax.plot_surface(range(x),range(y),self.dis[day,:,:],cmap = "coolwarm")
        return(self.dis[day,:,:])
