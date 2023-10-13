import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy.interpolate import interp1d

if __name__ == '__main__':

    
    apical_tension = np.arange(1.0, 2.0)
    dir_path =  os.path.dirname(os.path.realpath(__file__))
    A0 = [0.25, 0.5, 0.75, 1.0]
    radius = [2.5, 5.0]
    Ncells = [10, 10]
    exp_la = [2.5, 6.0]
    fig, ax = plt.subplots(len(radius), 1, figsize = (8, 8))
    lapical = np.zeros((len(radius), len(A0)))
    
    for r in range(len(radius)):
        A0_predict = np.arange(0, np.pi*radius[r]**2/Ncells[r] , 0.25)
        lapical_temp = np.zeros((len(A0), len(apical_tension)))
        lapical_formulas = (2*np.pi/Ncells[r])*np.sqrt(radius[r]**2 - Ncells[r]*A0_predict/np.pi)
        X_A0pre_int = interp1d(A0_predict, lapical_formulas, kind='cubic')
        X_A0_pre = np.linspace(min(A0_predict),max(A0_predict), 500)
        lapical_pred = X_A0pre_int(X_A0_pre)
        ax[r].plot(X_A0_pre, lapical_pred, color='purple', label='Theory')
        
        for area in range(len(A0)):
            for at in range(len(apical_tension)):
                file = open(f'{dir_path}/DifferentRadius/{radius[r]}/DifferentNcells/{Ncells[r]}/ApicalTension/{apical_tension[at]}/bundle_A0={A0[area]}.json', 'r')
                data = json.load(file)
                apical_y = data['Apical_y'][-1]
                apical_x = data['Apical_x'][-1]
                lapical_temp[area, at] = np.mean(data['Apical Lenght final']) 
        lapical[r, :] = np.mean(lapical_temp, axis=1)
        ax[r].scatter(A0, lapical[r,:], color='forestgreen', label = 'Simulation')
        ax[r].set_xlim(0, 2.0)
        #ax[r].vlines(x=1.4, ymax=exp_la[r] + 1, ymin=0, color='red')
        ax[r].hlines(y=exp_la[r], xmax=max(A0_predict), xmin=min(A0_predict), color='gray', linestyles='dashed')
        ax[r].set_title(f'Radius = {radius[r]}')
        ax[r].legend()

    fig.supylabel(f'Apical Length')
    fig.supxlabel(f'Preferred cell area (A0)')
    fig.savefig(f'{dir_path}/DifferentRadius/ALderivedZoom2.5-10.jpg')
