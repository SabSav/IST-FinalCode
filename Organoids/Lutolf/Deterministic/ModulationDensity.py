import matplotlib.pyplot as plt
import json
import numpy as np
import math
import os
import seaborn as sns
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


def fhc(x, Af, Ac, L, R, N, at):
    return ((np.pi**2)*Af/(2*L*Ac**2))*x**4 - 2*np.pi/Ac*(0.5 + R*np.pi*Af/(Ac*L))*x**3 + 2*np.pi/Ac*(R + (Af/L)*(0.5*N + np.pi*R**2/Ac))*x**2 - 2*np.pi*(N*Af*R/(Ac*L) + at)*x
if __name__ == '__main__':

    apical_tension = np.arange(1.0, 3.5, 0.5)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    l = 7.5
    rcup = 2.5
    Aflat = 0.5
    A_curved = 0.5
    Ncells = [20, 36, 44, 52, 60]
    density = [i/(2*np.pi*rcup + 2*l) for i in Ncells]
    hratioMean = np.zeros((len(apical_tension), len(Ncells)))
    hratioStd = np.zeros((len(apical_tension), len(Ncells)))
    mod = np.zeros((len(apical_tension), len(Ncells)))
    hc = np.zeros((len(apical_tension), len(Ncells)))
    hf = np.zeros((len(apical_tension), len(Ncells)))
    A_flat = 0.5
    A_curved = 0.5

    for at in range(len(apical_tension)):
        for ncell in range(len(Ncells)):
            file = open(f'{dir_path}/DifferentDensity/{Ncells[ncell]}/ApicalTension/{apical_tension[at]}/bundle.json')
            data = json.load(file)
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
            for vertex in range(Ncells[ncell]):
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
            for vertex in range(Ncells[ncell]):
                if (basal_x[vertex] == rcup and basal_y[vertex] not in y_flat_cross_plus) or (basal_x[vertex] == -rcup and basal_y[vertex] not in y_flat_cross_neg):
                    hflatSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))
                elif (basal_y[vertex] >= 0 and basal_x[vertex] not in x_curved_cross_plus) or (basal_y[vertex] <= -l and basal_x[vertex] not in x_curved_cross_neg):
                        hcurvedSimulation.append(math.sqrt((apical_x[vertex] - basal_x[vertex])**2 + (apical_y[vertex] - basal_y[vertex])**2))

            hratioMean[at, ncell] = np.mean(hflatSimulation)*np.mean([i**(-1) for i in hcurvedSimulation])
            hratioStd[at, ncell] = np.mean([i**2 for i in hflatSimulation])*np.mean([(i**(-1))**2 for i in hcurvedSimulation]) - (np.mean(hflatSimulation)**2)*(np.mean([i**(-1) for i in hcurvedSimulation]))**2
            mod[at, ncell] = -hratioMean[at, ncell] + 1
            res = minimize(fhc, x0=h, args=(A_flat, A_curved, l, rcup, Ncells[ncell], apical_tension[at]))
            hc[at, ncell] = res.x
            hf[at, ncell] = 0.5*Ncells[ncell]*A_flat/l - (A_flat*np.pi/(l*A_curved))*(rcup*hc[at, ncell] - 0.5*hc[at, ncell]**2)
    
    sns.set_style("whitegrid")
    color = []
    legends = []
    np.random.seed(0)
    for i in range(len(apical_tension)):
        color.append('#%06X' % np.random.randint(0, 0xFFFFFF))
        legends.append(mpatches.Patch(color=color[i], label=f'Γa = {apical_tension[i]}'))
    plt.figure(figsize=(8,6))
    for at in range(len(apical_tension)):
        X_r_hratio = interp1d(density, hf[at, :]/hc[at, :], kind='cubic')
        X_r = np.linspace(min(density), max(density), 500)
        hratio = X_r_hratio(X_r)
        plt.plot(X_r, -hratio+1, color=color[at])
        plt.errorbar(density, mod[at, :], yerr=hratioStd[at, :], fmt='o', color=color[at], label = f'Simulation')
    line = Line2D([0], [0], label='Theory', color='k')
    point = Line2D([0], [0], color='k', marker='o', markersize=5,ls="", label='Simulation')
    line_exp = Line2D([0], [0], label='Experiment', color='gray', linestyle='dashed')
    legends.append(line)
    legends.append(point)
    legends.append(line_exp)
    plt.hlines(y=0.6, xmin=min(density), xmax=1.6, linestyles='dashed', colors='gray')
    plt.hlines(y=0.7, xmin=min(density), xmax=1.6, linestyles='dashed', colors='gray')
    plt.legend(handles=legends, loc='upper right', framealpha=0.1)
    plt.xlabel('Density')
    plt.ylabel(r'$\Omega$', fontsize=15)
    plt.savefig(f'{dir_path}/DifferentDensity/ModDensity.jpg')



