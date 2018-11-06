# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 18:03:26 2018

@author: daniel.farber
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 12:03:59 2018

@author: daniel.farber
"""
#run model2d
#lateblight parameters from Andrade Piedra 2005
import model_2d_102818 as mlg #run the script containing the model (without calling the model)
import pickle #to save dictionaries
distribution = "inverse power" #dispesal kernel
parameters = [2.2, 25] #from lm(log10(Mean incidence) ~ log10(n+25),data=dg_df_tibble)
x = 101; y = 11; spacing = 1 #x-axis of field; y-axis of field; #row-spacing (e.g. "spacing = 1" means every row contains plants)
T = 40; l_period = 4; i_period = 16 #days of the epidemic; days of the latent period; days of the infectious period
dmfr = [1,2.5,5] #daily multiplication factor rate (the number of new infections created by each infections infection)
sites = 1000; #maximum number of infections per cell
initial_sev = [.01,.05,.1] #initial prevalence of latent infections (fraction > 0 and <=1)
'''think of focus as [row, column] not [x,y]'''
focus = [5,5] #[row(s),column(s)] of initial infections
lesion_growth = [0,1,2]
prod = [2.5,.05,1]

'''The Automated Way'''
alphabet = []
for letter in range(97,123):
    alphabet.append(chr(letter))
epi_str = 'Epidemic_217_5'
epi_list = []
pickle_list = []
from itertools import product
prod = list(product([0,1,2],repeat = 3))
for i in range(26):
     epi_list.append(epi_str + alphabet[i])
     pickle_list.append(epi_str + alphabet[i] +'.pickle')
for i in range(26):
     print(prod[i])
     print(epi_list[i])
     epi_list[i] =  mlg.model_2d(x,y,T,spacing,l_period,i_period,dmfr[prod[i][0]],sites,initial_sev[prod[i][1]],focus,distribution,parameters,lesion_growth[prod[i][2]])
     pickle_out = open(pickle_list[i],'wb')
     pickle.dump(epi_list[i],pickle_out)
     pickle_out.close()
     epi_list[i] = 0 


