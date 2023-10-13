import matplotlib.pyplot as plt
import json
import numpy as np
import os


if __name__ == '__main__':

    apical_tension = np.arange(1.0, 9.0)
    radius = [2.5, 5.0]
    A0=0.75
    dir_path = os.path.dirname(os.path.realpath(__file__))
    Ncells = 5
    fig, ax = plt.subplots(len(apical_tension), len(radius), figsize=(15, 40))
    fig.subplots_adjust(hspace=0, wspace=0)
    for at in range(len(apical_tension)):
        for rad in range(len(radius)):
            file = open(f'{dir_path}/DifferentRadius/{radius[rad]}/DifferentNcells/{Ncells}/ApicalTension/{apical_tension[at]}/bundle_A0={A0}.json')
            data = json.load(file)
            r = data['Rai']
            apical_y = data['Apical_y'][-1]
            apical_x = data['Apical_x'][-1]
            basal_x = data['Basal_x'][-1]
            basal_y = data['Basal_y'][-1]
            h = data['Thickness']
            x = []
            y = []
            for i in range(Ncells):
                x.append(apical_x[i])
                x.append(basal_x[i])
                y.append(apical_y[i])
                y.append(basal_y[i])
            for i in range(0, len(x),2):
                ax[at, rad].plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
            ax[at, rad].plot(radius[rad]*np.cos(np.arange(0, 2*np.pi, 1e-3)), radius[rad]*np.sin(np.arange(0, 2*np.pi, 1e-3)), color='darkred', label='basal')
            ax[at, rad].set_xlim(-5.5, 5.5)
            ax[at, rad].set_ylim(-5.5, 5.5)
            ax[at, rad].plot(apical_x[:], apical_y[:], color='navy', label='apical')
            ax[at, rad].plot([apical_x[0], apical_x[-1]], [apical_y[0], apical_y[-1]], color="navy")
            ax[len(apical_tension)-1, rad].set_title(f'{radius[rad]}', y=0, pad=-12)
            ax[at, 0].text(-7,0,f'{apical_tension[at]}', verticalalignment='center', rotation=-270)
            ax[at, rad].axis('off')
            
    fig.supxlabel('Radius',  y=np.sqrt(50) / 100 * 2.2 -7e-2, fontsize=18)
    fig.supylabel('Apical Tension', fontsize=18)
    fig.savefig(f'{dir_path}/DifferentRadius/FinalConfigs_A0={A0}.jpg', bbox_inches='tight')