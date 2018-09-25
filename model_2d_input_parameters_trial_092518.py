# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 14:56:32 2018

@author: daniel.farber
"""
#run model2d
import model2d #run the script containing the model (without calling the model)
distribution = "inverse power" #dispesal kernel
parameters = 2 #shape parameter
x = 51; y = 51; spacing = 1 #x-axis of field; y-axis of field; #row-spacing (e.g. "spacing = 1" means every row contains plants)
T = 75; l_period = 15; i_period = 15 #days of the epidemic; days of the latent period; days of the infectious period
dmfr = 5 #daily multiplication factor rate (the number of new infections created by each infections infection)
sites = 1000; #maximum number of infections per cell
initial_sev = .1 #initial prevalence of latent infections (fraction > 0 and <=1)
focus = [25,25] #[row(s),column(s)] of initial infections

epidemic = model2d.model_2d(x,y,T,spacing,l_period,i_period,dmfr,sites,initial_sev,focus,distribution,parameters)
