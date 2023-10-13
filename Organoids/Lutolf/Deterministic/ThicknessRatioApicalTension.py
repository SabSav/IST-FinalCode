import matplotlib.pyplot as plt
import json
import numpy as np
import math
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import os
import seaborn as sns

def fhc(x, Af, Ac, L, R, N, at):
    return ((np.pi**2)*Af/(2*L*Ac**2))*x**4 - 2*np.pi/Ac*(0.5 + R*np.pi*Af/(Ac*L))*x**3 + 2*np.pi/Ac*(R + (Af/L)*(0.5*N + np.pi*R**2/Ac))*x**2 - 2*np.pi*(N*Af*R/(Ac*L) + at)*x

if __name__ == '__main__':

    apical_tension = np.arange(1.0, 4.5, 0.5)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    hratioMean = []
    hratioStd = []
    Ncells = 20
    A_flat = 0.5
    A_curved = 0.5
    for at in apical_tension:
        file = open(f'{dir_path}/DifferentDensity/{Ncells}/ApicalTension/{at}/bundle.json')
        data = json.load(file)
        rcup = data['Rcups']
        l = data['Length']
        apical_y = data['Apical_y'][-1]
        apical_x = data['Apical_x'][-1]
        basal_x = data['Basal_x'][-1]
        basal_y = data['Basal_y'][-1]
        h = data['Thickness']
        hcurvedSimulation = []
        hflatSimulation = []
        hcross = []
        y_flat_plus = []
        y_flat_neg = []
        x_curved_plus = []
        x_curved_neg = []    
        for vertex in range(Ncells):
            if basal_x[vertex] == rcup and -l < basal_y[vertex] < 0 :
                y_flat_plus.append(basal_y[vertex])
            elif basal_x[vertex] == -rcup and -l < basal_y[vertex] < 0:
                y_flat_neg.append(basal_y[vertex])
            else:
                if basal_y[vertex] >= 0:
                    x_curved_plus.append(basal_x[vertex])
                elif basal_y[vertex] <= -l:
                    x_curved_neg.append(basal_x[vertex])
        y_flat_cross_plus = [max(y_flat_plus), min(y_flat_plus)]
        y_flat_cross_neg = [max(y_flat_neg), min(y_flat_neg)]
        x_curved_cross_neg = [max(x_curved_neg), min(x_curved_neg)]
        x_curved_cross_plus = [max(x_curved_plus), min(x_curved_plus)]
        for vertex in range(Ncells):
            if len(y_flat_neg) > 2 and len(y_flat_plus) > 2:
                if basal_y[vertex] in y_flat_cross_plus or basal_y[vertex] in y_flat_cross_neg:
                    hcross.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                elif basal_x[vertex] in x_curved_cross_plus or basal_x[vertex] in x_curved_cross_neg:
                    hcross.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                elif  (basal_x[vertex] == rcup and basal_y[vertex] not in y_flat_cross_plus) or (basal_x[vertex] == -rcup and basal_y[vertex] not in y_flat_cross_neg):
                    hflatSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                else:
                    hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
            else:
                if basal_x[vertex] == rcup or basal_x[vertex] == -rcup:
                    hflatSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                elif basal_x[vertex] in x_curved_cross_neg or basal_x[vertex] in x_curved_cross_plus:
                    basal_x[vertex] = basal_x[vertex]
                else:
                    hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
        
        hratioMean.append(np.mean(hflatSimulation)*np.mean([i**(-1) for i in hcurvedSimulation]))
        hratioStd.append(np.mean([i**2 for i in hflatSimulation])*np.mean([(i**(-1))**2 for i in hcurvedSimulation]) - (np.mean(hflatSimulation)**2)*(np.mean([i**(-1) for i in hcurvedSimulation]))**2)
       
    hc = np.zeros(len(apical_tension))
    hf = np.zeros(len(apical_tension))
    
    for i in range(len(apical_tension)):
        res = minimize(fhc, x0=h, args=(A_flat, A_curved, l, rcup, Ncells, apical_tension[i]))
        hc[i] = res.x
        hf[i] = 0.5*Ncells*A_flat/l - (A_flat*np.pi/(l*A_curved))*(rcup*hc[i] - 0.5*hc[i]**2)

    sns.set_style("whitegrid")
    X_at_hratio = interp1d(apical_tension, hf/hc, kind='linear')
    X_at = np.linspace(min(apical_tension), max(apical_tension), 500)
    hratio = X_at_hratio(X_at)
    plt.figure()
    plt.plot(X_at, hratio, color='purple', label='Theory')
    plt.errorbar(apical_tension, hratioMean, yerr=hratioStd, fmt='o', color = 'forestgreen', label = f'Simulation (Ncells = {Ncells})')
    plt.hlines(y=0.4, xmax=max(apical_tension), xmin=min(apical_tension), colors='gray', linestyles='dashed')
    plt.hlines(y=0.3, xmax=max(apical_tension), xmin=min(apical_tension), colors='gray', linestyles='dashed')
    plt.legend()
    plt.xlabel('Apical Tension')
    plt.ylabel(r'$\frac{h_{f}}{h_{c}}$', fontsize=15)
    plt.savefig(f'{dir_path}/DifferentDensity/{Ncells}/ApicalTension/ThicknessRatio.jpg')

    plt.figure()
    hratioMean_mod = [-i + 1 for i in hratioMean]
    plt.plot(X_at, -hratio+1, color='purple', label='Theory')
    plt.errorbar(apical_tension, hratioMean_mod, yerr=hratioStd, fmt='o', color = 'forestgreen', label = f'Simulation (Ncells = {Ncells})')
    plt.hlines(y=0.6, xmin=min(apical_tension), xmax=max(apical_tension), linestyles='dashed', colors='gray', label='Experiment')
    plt.hlines(y=0.7, xmin=min(apical_tension), xmax=max(apical_tension), linestyles='dashed', colors='gray')
    plt.legend(loc='lower right')
    plt.xlabel(r'$\Gamma_a$', fontsize=15)
    plt.ylabel(r'$\Omega$', fontsize=15)
    plt.savefig(f'{dir_path}/DifferentDensity/{Ncells}/ApicalTension/Modulation.jpg')




    
    