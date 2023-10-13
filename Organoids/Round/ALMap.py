import numpy as np
import matplotlib.pyplot as plt
import json
import os


if __name__ == '__main__':
    
    apical_tension = np.arange(1.0, 9.0)
    dir_path =  os.path.dirname(os.path.realpath(__file__))
    Ncells = 10
    A0 = [0.25, 0.5, 0.75, 1.0]
    radius = [10.0]
    fig, ax = plt.subplots(1, len(radius), figsize=(10, 4))
    for r in range(len(radius)):
        lapical = np.zeros((len(A0), len(apical_tension)))
        for area in range(len(A0)):
            for at in range(len(apical_tension)):
                thick = np.zeros(Ncells)
                file = open(f'{dir_path}/DifferentRadius/{radius[r]}/DifferentNcells/{Ncells}/ApicalTension/{apical_tension[at]}/bundle_A0={A0[area]}.json', 'r')
                data = json.load(file)
                h = data['Thickness']
                apical_y = data['Apical_y'][-1]
                apical_x = data['Apical_x'][-1]
                lapical[area, at] = np.mean(data['Apical Lenght final'])
                
        centers = [min(apical_tension), max(apical_tension), max(A0), min(A0)]
        dx = 1.0
        dy = 0.25
        extent = [centers[0]-dx/2, centers[1]+dx/2, centers[2]+dy/2, centers[3]-dy/2]
        plot = ax.imshow(lapical, interpolation='none', extent=extent, cmap=plt.cm.jet)  
        fig.colorbar(plot, ax=ax)
        ax.set_xticks(apical_tension)
        ax.set_yticks(np.flip(A0))
        ax.set_title(f'r={radius[r]}')
    fig.suptitle(f'Apical Length')
    fig.supxlabel(r'$\frac{\Gamma_a}{\Gamma_l}$', fontsize=15)
    fig.supylabel(r'$A_0$', fontsize=15)
    fig.savefig(f'{dir_path}/DifferentRadius/{radius[r]}/DifferentNcells/{Ncells}/ALMap.jpg')
