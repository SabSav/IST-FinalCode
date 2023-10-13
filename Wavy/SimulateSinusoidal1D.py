import Sinusoidal1D
import numpy as np
import matplotlib.pyplot as plt
import json
import os
from scipy.interpolate import interp1d


if __name__ == '__main__':

    dir_path =  os.path.dirname(os.path.realpath(__file__))

    np.random.seed(0)
    config = Sinusoidal1D.Energy(N_cells=30, alpha=0.201, beta=0.2, h=1, wavelength=1/3, apical_tension=0.6, compress_cell=20, lateral_tension=1, dt=1e-2)
    print(len(config.basal))
    plt.figure()
    x = []
    y = []
    for i in range(config.N_cells+1):
        x.append(config.apical[i,0])
        x.append(config.basal[i,0])
        y.append(config.apical[i,1])
        y.append(config.basal[i,1])
    for i in range(0, len(x),2):
        plt.plot(x[i:i+2], y[i:i+2], color="black", linewidth=0.5)
    X_Y_basal = interp1d(config.basal[:, 0], config.basal[:, 1], kind='cubic')
    Xbasal = np.linspace(np.min(config.basal[:, 0]), np.max(config.basal[:, 0]), 500)
    Y_basal = X_Y_basal(Xbasal)
    X_Y_apical = interp1d(config.apical[:, 0], config.apical[:, 1], kind='cubic')
    Xapical = np.linspace(np.min(config.apical[:, 0]), np.max(config.apical[:, 0]), 500)
    Y_apical = X_Y_apical(Xapical)
    plt.plot(Xbasal, Y_basal, color="black")
    plt.plot(Xapical, Y_apical, color="red")
    plt.xlim(0, config.N_cells)
    plt.ylim(-4, 6)
    plt.legend()
    plt.savefig(f'{dir_path}/SubstrateAmplitude/{config.beta}/InitialConfig.jpg')
    plt.close()


    config.areas = config.evaluate_areas() 
    config.A0[:] = 1
    plt.figure()
    bins = np.histogram_bin_edges(config.areas, 'doane')
    plt.hist(config.areas, bins=bins)
    plt.axvline(x = 1, color = 'red')
    plt.title(f'Initial areas distribution a = {config.apical_tension}, K = {config.compress_cell}')
    plt.savefig(f'{dir_path}/Wavelenght/{config.wavelength*config.N_cells}/InitialAreas.jpg')
    plt.close()

    config.energy()
    engy=[]
    index = 0
    apical_time_x = np.zeros((int(10**2/(config.dt)), config.N_cells+1))
    apical_time_y = np.zeros((int(10**2/(config.dt)), config.N_cells+1))
    basal_time_x = np.zeros((int(10**2/(config.dt)), config.N_cells+1))
    basal_time_y = np.zeros((int(10**2/(config.dt)), config.N_cells+1))
    for t in np.arange(0, 10**2, config.dt):
        engy.append(float(config.E))
        config.advance()
        index += 1
        print(f'a = {config.apical_tension}, K = {config.compress_cell} t = {t} , Energy = {config.E}')    

    bundle = config.export_dict()
    bundle['Wavelength'] = config.wavelength*config.N_cells
    bundle['ApicalTension'] = config.apical_tension
    bundle['Thickness'] = config.h
    bundle['BasalAmplitude'] = config.beta
    bundle['ApicalAmplitude'] = config.alpha
    bundle['Energy'] = engy
    bundle['Apical_x'] = config.apical[:, 0].tolist()
    bundle['Apical_y'] = config.apical[:, 1].tolist()
    bundle['Basal_x'] = config.basal[:, 0].tolist()
    bundle['Basal_y'] = config.basal[:, 1].tolist()

    with open(f'{dir_path}/SubstrateAmplitude/{config.beta}//Time.json', 'w') as file:
        json.dump(bundle, file, indent=12)    