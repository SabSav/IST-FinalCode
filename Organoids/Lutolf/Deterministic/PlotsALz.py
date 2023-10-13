import json
import numpy as np
import os
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


if __name__ == '__main__':

    l = 7.5
    Ncells=36
    dir_path = os.path.dirname(os.path.realpath(__file__))
    rcup = 2.5

    for apical_tension in np.arange(1.0, 8.5, 0.5):
        file = open(f'{dir_path}/{Ncells}/ApicalTension/{apical_tension}/bundle.json')
        z = []
        curvature = []
        data = json.load(file)
        h = data['Thickness']
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
        apical_y = data['Apical_y']
        apical_x = data['Apical_x']
        basal_x = data['Basal_x']
        basal_y = data['Basal_y']
        time = len(apical_y)
        lapical = np.zeros((time, len(z)))
        l_lat = np.zeros(len(z))
        for vertex in range(len(z)):
            l_lat[vertex] = math.sqrt((apical_x[-1][vertex] - basal_x[-1][vertex])**2 + (apical_y[-1][vertex] - basal_y[-1][vertex])**2)
            for t in range(time):
                lapical[t, vertex] = math.sqrt((apical_x[t][vertex] - apical_x[t][(vertex+1)%Ncells])**2 + (apical_y[t][vertex] - apical_y[t][(vertex+1)%Ncells])**2)
        for i in z:
            if 0 <= i < rcup:
                curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
            elif -l <= i <= 0:
                curvature.append(0)
            else:
                curvature.append(0.5*(1/(rcup - np.mean(l_lat)) + 1/(rcup)))
        z = np.array([0.5*(i - max(z))/(min(z) - max(z)) for i in z])        
        fig, ax = plt.subplots(2, 1, figsize=(20, 10))
        fig_time, ax_time = plt.subplots(2, 1, figsize=(20, 10))
        X_theta = np.linspace(z.min(), z.max(), 500)
        for t in range(time):
            X_Y_cubic_lapical = interp1d(z, lapical[t][0:len(z)], kind='cubic')
            Y_lapical = X_Y_cubic_lapical(X_theta)
            ax_time[1].plot(X_theta, Y_lapical, label=t)
        ax[1].plot(X_theta, Y_lapical, color='peru')
        ax[1].set_ylabel('Apical Length')
        ax_time[1].set_ylabel('Apical Length')    
        ax[0].step(z, curvature)
        ax[0].set_ylabel('Curvature')
        ax_time[0].step(z, curvature)
        ax_time[0].set_ylabel('Curvature')
        fig_time.legend()
        fig_time.supxlabel('z')
        fig_time.suptitle(f'Length flat side = {l}')
        fig_time.savefig(f'{dir_path}/{Ncells}/ApicalTension/{apical_tension}/LapicalTimeZcoord.jpg')
        fig.supxlabel('z')
        fig.suptitle(f'Length flat side = {l}')
        fig.savefig(f'{dir_path}/{Ncells}/ApicalTension/{apical_tension}/LapicalZcoord.jpg')
                
 
        