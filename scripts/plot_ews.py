import pickle

import pandas as pd
from numpy import random

from early_warning import get_largest_component, export_het_csv, calculate_indicators, plot_het_indicator
from funcs import assign_node_numbers

###### make a csv file for singlepop analysis using the earlywarning R package
# (10 nodes with the highest het
# 0erozygosity in the last step of each replica)
# read the pickle file with RGG (d-0.6) data
# with open(f'/home/lab2/PycharmProjects/fragmentation/RGG, cor_asymm_sig05_d06.pickle', 'rb') as file:
#     cor = pickle.load(file)
# print('finish')
# length, nets, het_dist, het_mean, fst_dist, nouse, nouse2, nouse3 = cor
# frag='cor'
# print(het_dist)
# export_het_csv(het_dist, frag)



###### make data for metapop analysis, get the het for the largest component
# with open('/home/lab2/PycharmProjects/fragmentation/RGG, cor_asymm_sig05_d06.pickle', 'rb') as file:
#     raw = pickle.load(file)
#
# length, nets, het_dist, het_mean, fst_dist, nouse, nouse2, nouse3 = raw
# het = assign_node_numbers(het_dist)
# components = get_largest_component(nets)
# component_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])
# component_data = component_data.sort_values(by=['replica', 'step', 'node_number'])
# component_data.reset_index().to_csv('cor_d0.6_component_sig05.csv', index=False)


#### calculate indicators for metapop (largest component)
# cor = pd.read_csv('/home/lab2/PycharmProjects/fragmentation/scripts/cor_d0.6_component_sig05.csv')
# indicators = calculate_indicators(cor)
# indicators.to_csv('indicators_metapop_sig05.csv', index=False)



###### plot het+indicators of metapopilation data-change y label
# random.seed(1)
# cor = pd.read_csv('/home/lab2/PycharmProjects/fragmentation/cor_d0.6_r100_component.csv')
# indicators = pd.read_csv('/home/lab2/PycharmProjects/fragmentation/indicators_metapop.csv')
# plot_het_indicator(cor, indicators, indicator='skew', n_samples=10)

###### plot het+indicators of single population data-change y label
# random.seed(1)
# cor = pd.read_csv('/home/lab2/PycharmProjects/fragmentation/cor_het.csv')
# indicators = pd.read_csv('/home/lab2/PycharmProjects/fragmentation/ASYMindicators_singlepop_25.csv')
# plot_het_indicator(cor, indicators, indicator='kurt', n_samples=10)