import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy.interpolate import interp1d


def evaluate_areas(Ncells, basal_x, basal_y, apical_x, apical_y):
    areas = np.zeros((Ncells, 1))
    j = 0
    for vertex in range(Ncells):
        x_coord = np.array([basal_x[vertex], basal_x[vertex+1], apical_x[vertex+1], apical_x[vertex]])
        y_coord = np.array([basal_y[vertex], basal_y[vertex+1], apical_y[vertex+1], apical_y[vertex]])
        A_cell = 0
        for i in range(len(x_coord)):
            A_cell += x_coord[i]*y_coord[(i+1) % len(x_coord)] - x_coord[(i+1) % len(x_coord)]*y_coord[i]
        areas[j] = np.abs(A_cell)*0.5
        j+=1

    return areas
    

if __name__ == '__main__':

    a = 2.0
    dir_path =  os.path.dirname(os.path.realpath(__file__))
    file = open(f'{dir_path}/ApicalTension/{a}/Time.json', 'r') 
    data = json.load(file)
    Ncells = data['Number of cells']
    h = data['Thickness']
    alpha = data['ApicalAmplitude']
    beta = data['BasalAmplitude']
    apical_y = data['Apical_y']
    apical_x = data['Apical_x']
    basal_x = data['Basal_x']
    basal_y = data['Basal_y']
    wavelength = data['Wavelength']
    energy = data['Energy']

    plt.figure()
    x = []
    y = []
    for i in range(Ncells+1):
        x.append(apical_x[i])
        x.append(basal_x[i])
        y.append(apical_y[i])
        y.append(basal_y[i])
    for i in range(0, len(x),2):
        plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
    X_Y_basal = interp1d(basal_x, basal_y, kind='cubic')
    Xbasal = np.linspace(min(basal_x), max(basal_x), 500)
    Y_basal = X_Y_basal(Xbasal)
    X_Y_apical = interp1d(apical_x, apical_y, kind='cubic')
    Xapical = np.linspace(min(apical_x), max(apical_x), 500)
    Y_apical = X_Y_apical(Xapical)
    plt.plot(np.arange(0, Ncells+1, 1e-3), alpha*np.cos(2*np.pi/(wavelength)*np.arange(0, Ncells + 1, 1e-3)) + h, color = 'gray',  linestyle='dashed', label='Initial Apical')
    plt.plot(Xbasal, Y_basal, color="black", label='Basal')
    plt.plot(Xapical, Y_apical, color="red", label= 'Apical')
    plt.xlim(0, Ncells)
    plt.ylim(-4, 6)
    plt.title(f'Apical Tension = {a}')
    plt.legend()
    plt.savefig(f'{dir_path}/ApicalTension/{a}/FinalConfig.jpg')
    plt.close()

    plt.figure()
    plt.plot(energy)
    plt.savefig(f'{dir_path}/ApicalTension/{a}/Energy.jpg')
    plt.close()
    

    plt.figure()
    a0 = 1
    areas = evaluate_areas(Ncells, basal_x, basal_y, apical_x, apical_y)
    bins = np.histogram_bin_edges(areas, 'doane')
    plt.hist(areas, bins=bins)
    plt.axvline(x = a0, color = 'red')
    plt.title(f'Final areas distribution')
    plt.savefig(f'{dir_path}/ApicalTension/{a}/FinalAreas.jpg')
    plt.close()