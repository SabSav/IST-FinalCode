import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import os
import json
import numpy as np
import math 
import seaborn as sns
import pandas as pd

if __name__ == '__main__':  
    
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    l = 20
    sigma = 0.1
    Nsamples = 100
    A0 = 0.5
    dt = 1e-3
    z = []
    rcup = 2.5
    file = open(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/bundle-0.json')
    data = json.load(file)
    Ncells = data['Number of cells']
    Nf = data['Number of cells on flat'] + 1
    Nc = int(0.5*(0.5*Ncells - Nf))
    theta = np.linspace(0, 0.5*np.pi, Nc)
    for i in range(len(theta)):
        z.append(np.cos(theta[i])*rcup)
    vertex = len(z) - 1
    for i in range(vertex + 1, vertex + 1 + Nf):
        z.append(z[i - 1] - l/Nf)
    theta = np.flip(theta)
    for i in range(1, len(theta)):
        z.append(-l - np.cos(theta[i])*rcup)
    std = np.zeros(Ncells)
    lapical = np.zeros((Nsamples, Ncells))
    l_lat = np.zeros((Nsamples, Ncells))

    for j in range(Nsamples):
        file = open(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/bundle-{j}.json')
        data = json.load(file)  
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        for vertex in range(Ncells):
            lapical[j, vertex] = math.sqrt((apical_x[vertex] - apical_x[(vertex+1)%Ncells])**2 + (apical_y[vertex] - apical_y[(vertex+1)%Ncells])**2)
            l_lat[j, vertex] = math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2)

    lapicalmean = np.mean(lapical, axis=0)
    sns.set_style("whitegrid")  # Setting style(Optional)
    std = np.std(lapical, axis=0) 
    l_lat = np.mean(l_lat, axis=0)
    l_lat = np.mean(l_lat)

    curvature = []
    for i in z:
        if 0 <= i < rcup:
            curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
        elif -l <= i <= 0:
            curvature.append(0)
        else:
            curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
    z = np.array([(i - max(z))/(min(z) - max(z)) for i in z])  
    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    X_theta = np.linspace(z.min(), z.max(), 500)
    X_Y_cubic_lapical = interp1d(z, lapicalmean[0:len(z)], kind='cubic')
    X_Y_cubic_std = interp1d(z, std[0:len(z)], kind='cubic')
    Y_lapical = X_Y_cubic_lapical(X_theta)
    Y_std = X_Y_cubic_std(X_theta)
    ax[1].plot(X_theta, Y_lapical, color='navy')
    ax[1].fill_between(X_theta, Y_lapical - Y_std, Y_lapical + Y_std, facecolor='navy', alpha=0.3)
    ax[1].set_ylabel(r'$\overline{l_{apical}}$', fontsize=15)
    ax[0].step(z, curvature, color='peru')
    ax[0].set_ylabel('Curvature')
    fig.supxlabel('Position (z/L)')
    fig.suptitle(f'Length flat side = {l}, σ = {sigma} (over 100 samples)')
    #fig.suptitle(f'Length flat side = {l}, Number of cells = {Ncells}, σ = {sigma} (over 100 samples)')
    fig.savefig(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/LapicalZcoord.jpg')


    fig, ax = plt.subplots(1, 2, figsize=(8, 8))
    ax[1].plot(X_theta, Y_lapical, color='navy')
    ax[1].fill_between(X_theta, Y_lapical - Y_std, Y_lapical + Y_std, facecolor='navy', alpha=0.3)
    ax[1].set_xlabel('Position (z/L)')
    sns.distplot(x = lapicalmean, kde = True , vertical='True', color = 'lightgray', kde_kws=dict(linewidth = 4 , color = 'navy'), ax=ax[0])
    fig.supylabel(r'$\overline{l_{apical}}}$' f'(over 100 samples)')
    fig.suptitle(f'Length flat side = {l}, σ = {sigma}')
    fig.savefig(f'{dir_path}/A0cases/{A0}/DifferentShape/DifferentLength/{l}/sigma/{sigma}/LapicalDistrib.jpg')