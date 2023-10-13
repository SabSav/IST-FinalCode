import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy.interpolate import interp1d


if __name__ == '__main__':
    
    r = 5.0
    apical_tension = np.arange(1.0, 9.0)
    A0 = 0.75
    Ncells=20
    dir_path =  os.path.dirname(os.path.realpath(__file__))
    for at in apical_tension:
        file = open(f'{dir_path}/DifferentRadius/{r}/DifferentNcells/{Ncells}/ApicalTension/{at}/bundle_A0={A0}.json', 'r') 
        data = json.load(file)
        Ncells = data['Number of cells']
        h = data['Thickness']

        #plot initial configuration
        apical_y = data['Apical_y'][0]
        apical_x = data['Apical_x'][0]
        basal_x = data['Basal_x'][0]
        basal_y = data['Basal_y'][0]

        plt.figure()
        x = []
        y = []
        for i in range(Ncells):
            x.append(apical_x[i])
            x.append(basal_x[i])
            y.append(apical_y[i])
            y.append(basal_y[i])
            
        for i in range(0, len(x),2):
            plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)

        plt.plot(r*np.cos(2*np.pi*np.arange(0, Ncells, 1e-3)), r*np.sin(2*np.pi*np.arange(0, Ncells, 1e-3)), color="darkred")
        plt.plot((r-h)*np.cos(2*np.pi*np.arange(0, Ncells, 1e-3)), (r-h)*np.sin(2*np.pi*np.arange(0, Ncells, 1e-3)), color="navy")
        plt.scatter(basal_x , basal_y, color="darkred", label='Basal')
        plt.scatter(apical_x , apical_y, color="navy", label='Apical')
        plt.legend()
        plt.savefig(f'{dir_path}/DifferentRadius/{r}/DifferentNcells/{Ncells}/ApicalTension/{at}/Initial-{A0}.jpg')
        plt.close()

        #plot energy
        plt.figure()
        energy = data['Energy']
        plt.plot(energy)
        plt.savefig(f'{dir_path}/DifferentRadius/{r}/DifferentNcells/{Ncells}/ApicalTension/{at}/Energy-{A0}.jpg')
        plt.close()

        #plot final configuration and final cell areas distribution
        
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        plt.figure()
        x = []
        y = []
        for i in range(Ncells):
            x.append(apical_x[i])
            x.append(basal_x[i])
            y.append(apical_y[i])
            y.append(basal_y[i])
        for i in range(0, len(x),2):
            plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)

        plt.plot(r*np.cos(2*np.pi*np.arange(0, Ncells, 1e-3)), r*np.sin(2*np.pi*np.arange(0, Ncells, 1e-3)), color="darkred")
        plt.plot((r-h)*np.cos(2*np.pi*np.arange(0, Ncells, 1e-3)), (r-h)*np.sin(2*np.pi*np.arange(0, Ncells, 1e-3)), color="navy")
        plt.scatter(basal_x , basal_y, color="darkred", label='Basal')
        plt.scatter(apical_x , apical_y, color="navy", label='Apical')
        plt.title(f'Apical Tension = {at}')
        plt.legend()
        plt.savefig(f'{dir_path}/DifferentRadius/{r}/DifferentNcells/{Ncells}/ApicalTension/{at}/Final-{A0}.jpg')
        plt.close()