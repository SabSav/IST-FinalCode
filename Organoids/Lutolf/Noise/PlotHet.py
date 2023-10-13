import json
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from scipy.interpolate import interp1d

if __name__ == '__main__':

    l = 7.5
    sigmas = [0, 0.01, 0.05, 0.1, 0.2]
    Nsample = 100
    Ncells = 20
    dir_path = os.path.dirname(os.path.realpath(__file__))
    rcup = 2.5

    #deterministic
    file = open(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{0}/bundle.json')
    data = json.load(file)
    Nf = data['Number of cells on flat'] + 1
    Nc = int(0.5*(0.5*Ncells - Nf))
    h = data['Thickness']
    apical_y_0 = data['Apical_y'][-1]
    apical_x_0 = data['Apical_x'][-1]
    basal_x_0 = data['Basal_x'][-1]
    basal_y_0 = data['Basal_y'][-1]
    het_sigmas = np.zeros((len(sigmas), Ncells))
    yap = np.zeros(Ncells)
    het = np.zeros(Ncells)
    l_lat = np.zeros(Ncells)
    for vertex in range(Ncells):
        yap[vertex] = math.sqrt((apical_x_0[vertex] - apical_x_0[(vertex + 1)%Ncells])**2 + (apical_y_0[vertex] - apical_y_0[(vertex + 1)%Ncells])**2)
        l_lat[vertex] = math.sqrt((apical_x_0[vertex] - basal_x_0[vertex])**2 + (apical_y_0[vertex] - basal_y_0[vertex])**2) 
    for vertex in range(Ncells):
        het_sigmas[0, vertex] = 0.5*((yap[vertex] - yap[(vertex+1)%Ncells])**2 + (yap[vertex] - yap[(vertex-1)%Ncells])**2)
    l_lat = np.mean(l_lat)

    #noise
    for sig in range(1,len(sigmas)):
        yap = np.zeros((Nsample, Ncells))
        het = np.zeros((Nsample, Ncells))
        for j in range(Nsample):
            file = open(f'{dir_path}/DifferentDensity/{Ncells}/sigma/{sigmas[sig]}/bundle-{j}.json')
            data = json.load(file)
            apical_y = data['Apical_y'][-1]
            apical_x = data['Apical_x'][-1]
            for vertex in range(Ncells):
                yap[j, vertex] = math.sqrt((apical_x[vertex] - apical_x[(vertex+1)%Ncells])**2 + (apical_y[vertex] - apical_y[(vertex+1)%Ncells])**2)
            for vertex in range(Ncells):    
                het[j, vertex] = 0.5*((yap[j, vertex] - yap[j, (vertex+1)%Ncells])**2 + (yap[j, vertex] - yap[j, (vertex-1)%Ncells])**2)
                
        het_sigmas[sig, :] = np.mean(het, axis=0)

    het_noise = np.zeros((len(sigmas)-1, Ncells))
    for sig in range(1, len(sigmas)):
        het_noise[sig-1, :] = het_sigmas[sig, :] - het_sigmas[0, :]
    
    z = []
    theta = np.linspace(0, 0.5*np.pi, Nc)
    for i in range(len(theta)):
        z.append(np.cos(theta[i])*rcup)
    vertex = len(z) - 1
    for i in range(vertex + 1, vertex + 1 + Nf):
        z.append(z[i - 1] - l/Nf)
    theta = np.flip(theta)
    for i in range(1, len(theta)):
        z.append(-l - np.cos(theta[i])*rcup)
    curvature = []
    for i in z:
        if 0 <= i < rcup:
            curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
        elif -l <= i <= 0:
            curvature.append(0)
        else:
            curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
    z = np.array([(i - max(z))/(min(z) - max(z)) for i in z])  

    fig, ax = plt.subplots(len(sigmas), 1, figsize = (12, 12))
    ax[0].step(z, curvature, color='peru')
    ax[0].set_ylabel('Curvature')
    X_theta = np.linspace(z.min(), z.max(), 500)
    X_Y_cubic_het_det = interp1d(z, het_sigmas[0, 0:len(z)], kind='cubic')
    Y_het_det = X_Y_cubic_het_det(X_theta)
    for i in range(1, len(sigmas)):
        X_Y_cubic_het = interp1d(z, het_sigmas[i, 0:len(z)], kind='cubic')
        Y_het = X_Y_cubic_het(X_theta)
        X_Y_cubic_het_noise = interp1d(z, het_noise[i-1, 0:len(z)], kind='cubic')
        Y_het_noise = X_Y_cubic_het_noise(X_theta)
        ax[i].plot(X_theta, Y_het_det, color='darkblue', label='deterministic')
        ax[i].plot(X_theta, Y_het_noise, color='magenta', label='noise')
        ax[i].plot(X_theta, Y_het, color='red', label='stochastic')
        ax[i].legend()
        ax[i].set_ylabel('Heterogeneity')
        ax[i].set_title(f'σ = {sigmas[i]}')
    fig.supxlabel('Position (z/L)')
    fig.suptitle(f'Averages over {Nsample} organoids')
    fig.savefig(f'{dir_path}/DifferentDensity/{Ncells}/Het.jpg',bbox_inches='tight')